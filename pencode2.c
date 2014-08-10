/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Description:                                                                * 
 *                                                                             *
 * pencode2.c for Assignment 3 of COMP9319, 14s1                               *
 * A position-couscious pwt encoder that expects 2 command line arguments:     *
 *    argv[1]-the original file to read, argv[2]-the pwt file to write to.     *
 * This program encodes a pwt file in linear time, but is different from       *
 * pencode.c in that it uses temporaty external files to store suffix array    *
 * buckets and requires less memory (only the original input is in memory).    *
 *    1. Create 128 temporary files as buckets.                                *
 *    2. Extract all '[' positions in order, fill relevant L positions into    *
 *       corresponding temporary files. Write pwt chars corresponding to '['   *
 *       positions in the suffix array to output file.                         *
 *    3. For each char in alphabet, extract its S positions, sort them.        *
 *       combine them with L positions from temporary files, and fill relevant *
 *       L positions into larger buckets in temporary files. Write pwt chars   *
 *       corresponding to this bucket of the suffix array to output file.      *
 *                                                                             *
 * Written by Lei Qian, 3447461                                                *
 * Other source files, if any, one per line, starting on the next line:        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#define ASCII 256
unsigned char *original;
int file_size;
int nb_of_buckets;
int nb_of_new_buckets;
unsigned char *is_bucket_end;
int radix;
int compare(const void *a, const void *b);
int compare_char(int p1, int p2);
bool is_Small(int p);
int compare_position_after_shift(const void *p1, const void *p2);
void cumulate(int *occ, int *cumulate_occ);
void MSD_radixsort(int *positions, int start, int end);
void remaining_quicksort(int *positions, int start, int end);

int main(int argc, char **argv) {
    FILE *input = fopen(argv[1], "rb");
    fseek(input, 0, SEEK_END);
    file_size = ftell(input);
    int total_occ[ASCII] = {0};
    /* Read the whole input file into memory. */
    original = malloc(file_size);
    fseek(input, 0, SEEK_SET);
    fread(original, 1, file_size, input);
    fclose(input);
    /* Occurrence of each S-type char. */
    int small_occ[ASCII] = {0};
    FILE *des = fopen(argv[2], "wb");
    /* Count total_occ and small_occ */
    for (int i = 0; i < file_size; ++i) {
        total_occ[original[i]]++;
        if (is_Small(i))
            small_occ[original[i]]++;
    }
    /* File names and FILE pointers to temporary files as buckets. */
    char tmp_name[ASCII][10];
    FILE *tmp_bucket[ASCII];
    for (int c = 0; c < ASCII; ++c)
        if (total_occ[c] && c != '[') {
            sprintf(tmp_name[c], "%d.tmp", c);
            tmp_bucket[c] = fopen(tmp_name[c], "wb");
        }
    /* char_buffer stores pwt chars. */
    unsigned char *char_buffer = malloc(total_occ['[']);
    /* bucket_buffer stores positions in the suffix array, one bucket at a time. */
    int *bucket_buffer = malloc(sizeof(int) * total_occ['[']);
    /* bucket_index_S traces the indexes of each S position in a bucket. */
    int bucket_index_S = 0;
    /* Scan and read '[' positions into bucket '[' of suffix array. */
    for (int i = 0; i < file_size; ++i)
        if (original[i] == '[')
            bucket_buffer[bucket_index_S++] = i;
    /* bucket_index_L[] maintains the indexes of each L position in a bucket. */
    int bucket_index_L[ASCII] = {0};
    /* Use positions in '[' bucket to update L positions in following buckets, */
    /* then write chars corresponding to the '[' bucket of suffix array to output file.*/
    for (int i = 0; i < total_occ['[']; ++i) {
        int previous_position = (bucket_buffer[i] + file_size - 1) % file_size;
        if (is_Small(previous_position))
            continue;
        unsigned char bucket_number = original[previous_position];
        if (bucket_number == '[')
            continue;
        bucket_index_L[bucket_number]++;
        /* Write L positions derived from '[' bucket to their corresponding bucket files. */
        fwrite(&previous_position, sizeof(int), 1, tmp_bucket[bucket_number]);
        /* Write pwt chars based on positions in the '[' bucket of suffix array. */
        char_buffer[i] = original[(bucket_buffer[i] + file_size - 1) % file_size];
    }
    fwrite(char_buffer, 1, total_occ['['], des);
    free(bucket_buffer);
    free(char_buffer);


    /* For each char in alphabet: */
    /* (1) read L positions from temporary files; */
    /* (2) extract S positions and sort them; */
    /* (3) combine S positions with L positions, derive relevant L positions and fill them */
    /*     in following buckets in temporary files, or update current bucket in memory. */
    /* (4) write pwt chars corresponding to this bucket of the suffix array to output file. */
    for (int c = 0; c < ASCII; ++c) {
        /* Skip '[' bucket and empty buckets. */
        if (c == '[' || total_occ[c] == 0)
            continue;
        /* Close the bucket file (in write-mode) and reopen it in read-mode. */
        /* Then read L positions into suffix array. */
        fclose(tmp_bucket[c]);
        tmp_bucket[c] = fopen(tmp_name[c], "rb");
        bucket_buffer = malloc(sizeof(int) * total_occ[c]);
        char_buffer = malloc(total_occ[c]);
        fread(bucket_buffer, sizeof(int), total_occ[c], tmp_bucket[c]);
        /* Extract S positions with this char c and place them in suffix array, */
        /* following the existent L positions. */
        bucket_index_S = 0;
        for (int i = 0; i < file_size; ++i)
            if (original[i] == c && is_Small(i))
                bucket_buffer[total_occ[c] - small_occ[c] + (bucket_index_S++)] = i;
        
        /* Within each bucket except '[' bucket, radix sort all S positions. When radix sort */
        /* becomes inefficient, use quicksort to resolve remaining sub-buckets. */
        
        /* MSD radix sort each bucket from radix 1, i.e., comparing the char that is */
        /* one "digit" after the current position. */
        /* All elements in a bucket have identical beginning char, so radix sort should */
        /* begin with radix shift of at least 1. */
        nb_of_buckets = 1;
        radix = 1;
        /* The number of new buckets generated in the last round of radix sort. */
        /* Initialized as an arbitrarily large number to avoid triggering cut-off. */
        nb_of_new_buckets = 10000;
        /* Label the end of current alphabet bucket if it is not empty. */
        is_bucket_end = calloc(1, (total_occ[c] - 1) / 8 + 1);
        is_bucket_end[(total_occ[c] - 1) / 8] |= 1 << ((total_occ[c] - 1) % 8);
        /* Loop until all elements are resolved. */
        while (nb_of_buckets < small_occ[c]) {
            /* When efficiency of radix sort falls under threshold based on cut-off rule, use quick sort to resolve remaining sub-buckets. */
            if (1.0 * nb_of_new_buckets / small_occ[c] < 0.1 && nb_of_buckets > small_occ[c] * 0.7) {
                remaining_quicksort(bucket_buffer, total_occ[c] - small_occ[c], total_occ[c] - 1);
                nb_of_buckets = small_occ[c];
            }
            /* MSD Radix sort at the current radix. */
            else {
                nb_of_new_buckets = 0;
                MSD_radixsort(bucket_buffer, total_occ[c] - small_occ[c], total_occ[c] - 1);
                nb_of_buckets += nb_of_new_buckets;
                radix++;
                if (radix == file_size)
                    break;
            }
        }
        /* Write L positions derived from '[' bucket to their corresponding bucket files, */
        /* or update suffix array in memory. */
        for (int i = 0; i < total_occ[c]; ++i) {
            int previous_position = (bucket_buffer[i] + file_size - 1) % file_size;
            /* Skip S positions, which are already placed in the right location. */
            if (is_Small(previous_position))
                continue;
            unsigned char bucket_number = original[previous_position];
            int new_location = bucket_index_L[bucket_number]++;
            /* For L positions in following buckets, write to temporary files. */
            if (bucket_number > c && bucket_number != '[')
                fwrite(&previous_position, sizeof(int), 1, tmp_bucket[bucket_number]);
            /* For L positions in the current bucket, update suffix array in memory. */
            else if (bucket_number == c)
                bucket_buffer[new_location] = previous_position;
        }
        /* Write pwt chars based on positions in the current bucket of suffix array. */
        for (int i = 0; i < total_occ[c]; ++i)
            char_buffer[i] = original[(bucket_buffer[i] + file_size - 1) % file_size];
        fwrite(char_buffer, 1, total_occ[c], des);
        free(bucket_buffer);
        free(char_buffer);
        fclose(tmp_bucket[c]);
        remove(tmp_name[c]);
        free(is_bucket_end);
    }
    
    free(original);
    fclose(des);
    return EXIT_SUCCESS;
}

/* Resolve two positions based on the given radix. */
int compare(const void *a, const void *b) {
    for (int i = 0; i < file_size - radix; ++i) {
        int p1 = (*(int *)a + radix + i) % file_size;
        int p2 = (*(int *)b + radix + i) % file_size;
        if (original[p1] == '[') {
            if (original[p2] != '[')
                return -1;
            if (p1 > p2)
                return 1;
            else
                return -1;
        }
        if (original[p2] == '[')
            return 1;
        if (original[p1] > original[p2])
            return 1;
        if (original[p1] < original[p2])
            return -1;
    }
    return 1;
}

/* Compare two positions based on the chars on these positions. */
/* No further compare on following chars. */
/* Used to differentiate S and L positions. */
int compare_char(int p1, int p2) {
    if (original[p1] == '[') {
        if (original[p2] != '[')
            return -1;
        return p1 > p2 ? 1 : -1;
    }
    if (original[p2] == '[')
        return 1;
    if (original[p1] > original[p2])
        return 1;
    if (original[p1] < original[p2])
        return -1;
    return 0;
}

/* Test if a position is S-type. */
bool is_Small(int p) {
    for (int i = 0; i < file_size - 1; ++i) {
        if (compare_char(p + i, (p + i + 1) % file_size) > 0)
            return false;
        else if (compare_char(p + i, (p + i + 1) % file_size) < 0)
            return true;
    }
    return true;
}

/* Compare two '[' positions based on radix shift, used to quicksort '[' positions. */
int compare_position_after_shift(const void *p1, const void *p2) {
    return (*(int *)p1 + radix) % file_size > (*(int *)p2 + radix) % file_size ? 1 : -1;
}

/* Calculate the number of chars that are smaller than the current char. */
void cumulate(int *occ, int *result) {
    result[0] = occ['['];
    for (int i = 1; i < ASCII; ++i) {
        if (i == '[')
            result[i] = 0;
        else if (i == '[' + 1)
            result[i] = result[i - 2] + occ[i - 2];
        else
            result[i] = result[i - 1] + occ[i - 1];
    }
}

/* Radix sort part of an array of positions based on the given radix. */
void MSD_radixsort(int *positions, int start, int end) {
    int sub_bucket_size = 0;
    int occ[ASCII] = {0};
    /* Scan the array and find every sub-bucket that is not fully resolved. */
    /* At the same time, keep a record of occurrence of every char, which will */
    /* be used to determine their location when sorted in-place. */
    for (int i = start; i <= end; ++i) {
        /* A bucket end is encountered. */
        if ((is_bucket_end[i / 8] & 1 << (i % 8))) {
            /* Skip sub-buckets that are fully resolved. */
            if (sub_bucket_size == 0)
                continue;
            sub_bucket_size++;
            /* Quicksort sufficiently small sub-buckets to improve efficiency. */
            if (sub_bucket_size <= 1024 * 8) {
                nb_of_new_buckets += (sub_bucket_size - 1);
                for (int j = i - sub_bucket_size + 1; j < i; ++j)
                    is_bucket_end[j / 8] |= 1 << (j % 8);
                qsort(&positions[i - sub_bucket_size + 1], sub_bucket_size, sizeof(int), compare);
                for (int j = 0; j < ASCII; ++j)
                    occ[j] = 0;
                sub_bucket_size = 0;
                continue;
            }
            /* In-place bucket sort the sub-bucket. */
            int sorted_position;
            bool *sorted = calloc(1, sub_bucket_size);
            occ[original[(positions[i] + radix) % file_size]]++;
            int cumulate_occ[ASCII];
            int bucket_index[ASCII] = {0};
            cumulate(occ, cumulate_occ);
            for (int j = i - sub_bucket_size + 1; j <= i; ++j) {
                if (sorted[j - i + sub_bucket_size - 1])
                    continue;
                int pending_position = positions[j];
                unsigned char radix_char = original[(pending_position + radix) % file_size];
                sorted_position = i - sub_bucket_size + 1 + cumulate_occ[radix_char] + (bucket_index[radix_char]++);
                while (sorted_position > j) {
                    int tmp = positions[sorted_position];
                    sorted[sorted_position - i + sub_bucket_size - 1] = true;
                    positions[sorted_position] = pending_position;
                    /* Label the end of new sub-buckets. */
                    if (bucket_index[radix_char] == occ[radix_char] && !(is_bucket_end[sorted_position / 8] & 1 << (sorted_position % 8))) {
                        nb_of_new_buckets++;
                        is_bucket_end[sorted_position / 8] |= 1 << (sorted_position % 8);
                    }
                    pending_position = tmp;
                    radix_char = original[(pending_position + radix) % file_size];
                    sorted_position = i - sub_bucket_size + 1 + cumulate_occ[radix_char] + (bucket_index[radix_char]++);
                }
                sorted[sorted_position - i + sub_bucket_size - 1] = true;
                positions[sorted_position] = pending_position;
                if (bucket_index[radix_char] == occ[radix_char] && !(is_bucket_end[sorted_position / 8] & 1 << (sorted_position % 8))) {
                    nb_of_new_buckets++;
                    is_bucket_end[sorted_position / 8] |= 1 << (sorted_position % 8);
                }
            }
            /* Resolve the '[' sub-bucket using quicksort. */
            if (occ['['] > 1) {
                qsort(&positions[i - sub_bucket_size + 1], occ['['], sizeof(int), compare_position_after_shift);
                for (int j = i - sub_bucket_size + 1; j <= i - sub_bucket_size + occ['[']; ++j)
                    is_bucket_end[j / 8] |= 1 << (j % 8);
                nb_of_new_buckets += (occ['['] - 1);
            }
            free(sorted);
            for (int j = 0; j < ASCII; ++j)
                occ[j] = 0;
            sub_bucket_size = 0;
        }
        /* A bucket end is not yet encountered. */
        else {
            sub_bucket_size++;
            occ[original[(positions[i] + radix) % file_size]]++;
        }
    }
}

/* Quick sort part of an array of positions based on the given radix. */
void remaining_quicksort(int *positions, int start, int end) {
    int sub_bucket_size = 0;
    for (int i = start; i <= end; ++i) {
        if ((is_bucket_end[i / 8] & 1 << (i % 8))) {
            if (sub_bucket_size == 0)
                continue;
            sub_bucket_size++;
            qsort(&positions[i - sub_bucket_size + 1], sub_bucket_size, sizeof(int), compare);
            sub_bucket_size = 0;
        }
        else
            sub_bucket_size++;
    }
}
