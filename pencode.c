/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Description:                                                                * 
 *                                                                             *
 * pencode.c for Assignment 3 of COMP9319, 14s1                                *
 * A position-couscious pwt encoder that expects 2 command line arguments:     *
 *    argv[1]-the original file to read, argv[2]-the pwt file to write to.     *
 * This program encodes a pwt file in linear time as follows:                  *
 *    1. Extract all Small-type positions.                                     *
 *    2. Construct a suffix array and divide it into 128 buckets.              *
 *    3. Place all positions with char '[' into the first bucket based on the  *
 *       position values. These bucket is finalized.                           *
 *    4. Place the other S positions into corresponding buckets from the back  *
 *       ends.                                                                 *
 *    5. Radix sort S positions in each bucket. When appropriate, apply        *
 *       quicksort to speed up the resolving process. What is different from   *
 *       bencode.c include: (1) '[' is deemed smaller than any other char;     *
 *       (2) two '['s are resolved by their position.                          *
 *    6. Place L positions and fill the gaps in the suffix array.              * 
 *    7. Construct bwt file based on suffix array.                             *
 *                                                                             *
 * Written by Lei Qian, 3447461                                                *
 * Other source files, if any, one per line, starting on the next line:        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#define ASCII 256
unsigned char *original;
int file_size;
int nb_of_buckets;
int nb_of_new_buckets;
unsigned char *is_bucket_end;
int radix;
int compare(const void *a, const void *b);
int compare_char(int p1, int p2);
int compare_position_after_shift(const void *p1, const void *p2);
void cumulate(int *occ, int *cumulate_occ);
void inverse_cumulate(int *occ, int *inverse_cumulate_occ);
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
    /* Scan the input and label whether a char is S-type or L-type. */
    /* Default as not small. */
    /* Use a bit array to save space. */
    unsigned char *is_small = calloc(1, (file_size - 1) / 8 + 1);
    int small_count = 0;
    /* Total occurrence of each char. */
    total_occ[original[0]]++;
    /* Occurrence of each S-type char. */
    int small_occ[ASCII] = {0};
    int last_pending_position = 0;
    for (int i = 1; i < file_size; ++i) {
        total_occ[original[i]]++;
        if (compare_char(i - 1, i) < 0) {
            for (int j = last_pending_position; j < i; ++j) {
                is_small[j / 8] |= 1 << (j % 8);
                small_count++;
                small_occ[original[j]]++;
            }
            last_pending_position = i;
        }
        else if (compare_char(i - 1, i) > 0)
            last_pending_position = i;
    }
    /* Since the last char is no special $, we need to compare it with the first char to determine if it is S or L. */
    if (compare_char(file_size - 1, 0) < 0 || (compare_char(file_size - 1, 0) == 0 && is_small[0] & 1)) {
        for (int j = last_pending_position; j < file_size; ++j) {
            is_small[j / 8] |= 1 << (j % 8);
            small_count++;
            small_occ[original[j]]++;
        }
    }

    int *suffix_array = malloc(sizeof(int) * file_size);
    for (int i = 0; i < file_size; ++i) {
        suffix_array[i] = -1;
    }
    /* Build the C[] array. */
    int cumulate_total_occ[ASCII] = {0};
    cumulate(total_occ, cumulate_total_occ);
    /* Bucket sort all S substring positions for the first round. */
    int bucket_index[ASCII] = {0};
    /* If the last char is '[', it must be placed at the end of bucket '[', */
    /* no matter it is S-type or L-type. */
    if (original[file_size - 1] == '[')
        suffix_array[cumulate_total_occ['['] + total_occ['['] - 1 - (bucket_index['[']++)] = file_size - 1;
    /* Else, place the last char in the right location. */
    else if ((is_small[(file_size - 1) / 8] & 1 << ((file_size - 1) % 8)))
        suffix_array[cumulate_total_occ[original[file_size - 1]] + total_occ[original[file_size - 1]] - 1 - (bucket_index[original[file_size - 1]]++)] = file_size - 1;
    /* Bucket sort all S substring positions for the first round. */
    /* Scan the input backward because we want '[' positions to be in acending order. */
    for (int i = file_size - 2; i >= 0; --i) {
        /* An L position encountered, skip. */
        if ((is_small[i / 8] & 1 << (i % 8)) == 0) {
            continue;
        }
        /* Place an S position into its bucket through the back end. */
        /* [-1, -1, ...., S2, S1, S0]*/
        suffix_array[cumulate_total_occ[original[i]] + total_occ[original[i]] - 1 - (bucket_index[original[i]]++)] = i;
    }

    /* Within each bucket except '[' bucket, radix sort all S positions. When radix sort */
    /* becomes inefficient, use quicksort to resolve remaining sub-buckets. */

    /* Label the end of sub-bucket using a bit array. */
    is_bucket_end = calloc(1, (file_size - 1) / 8 + 1);
    for (int i = 0; i < ASCII; ++i) {
        /* '[' bucket is already sorted by position, so skip it. */
        if (i == '[' || small_occ[i] == 0)
            continue;
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
        if (total_occ[i] > 0) {
            int tmp = cumulate_total_occ[i] + total_occ[i] - 1;
            is_bucket_end[tmp / 8] |= 1 << (tmp % 8);
        }
        /* Loop until all elements are resolved. */
        while (nb_of_buckets < small_occ[i]) {
            /* When radix sort becomes inefficient, quicksort the remaining sub-buckets. */
            if (1.0 * nb_of_new_buckets / small_occ[i] < 0.1 && nb_of_buckets > small_occ[i] * 0.7) {
                /* When efficiency of radix sort falls under threshold based on cut-off rule, use quick sort to resolve remaining sub-buckets. */
                remaining_quicksort(suffix_array, cumulate_total_occ[i] + total_occ[i] - small_occ[i], cumulate_total_occ[i] + total_occ[i] - 1);
                nb_of_buckets = small_occ[i];
            }
            /* MSD Radix sort at the current radix. */
            else {
                nb_of_new_buckets = 0;
                MSD_radixsort(suffix_array, cumulate_total_occ[i] + total_occ[i] - small_occ[i], cumulate_total_occ[i] + total_occ[i] - 1);
                nb_of_buckets += nb_of_new_buckets;
                radix++;
                if (radix == file_size)
                    break;
            }
        }
    }
    for (int i = 0; i < ASCII; ++i) {
        bucket_index[i] = 0;
    }
    /* In case all chars are identical so there are no S positions, */
    /* pick 0 as the first number of suffix array. */
    if (suffix_array[0] == -1)
        suffix_array[0] = 0;
    /* Place all L positions into suffix array from the first to the last bucket. */
    for (int i = 0; i < file_size; ++i) {
        int previous_position = suffix_array[i] ? suffix_array[i] - 1 : file_size - 1;
        /* Skip S positions, which are already placed in the right location. */
        if (is_small[previous_position / 8] & 1 << (previous_position % 8))
            continue;
        unsigned char bucket_number = original[previous_position];
        /* Place an L position into its bucket through the front end. */
        /* [L0, L6, L7, ...., S, S, S] */
        /* [L1, L3, L8, ...., S, S, S] */
        /* [L2, L4, L5, ...., S, S, S] */
        suffix_array[cumulate_total_occ[bucket_number] + (bucket_index[bucket_number]++)] = previous_position;
    }
    
    /* Output bwt file. */    
    int buffer_size = 1024 * 512;
    unsigned char *char_buffer = malloc(buffer_size);
    int output_times = (file_size - 1) / buffer_size + 1;
    FILE *destination = fopen(argv[2], "wb");
    for (int i = 0; i < output_times; ++i) {
        if (i < output_times - 1) {
            for (int j = 0; j < buffer_size; ++j)
                char_buffer[j] = original[(suffix_array[i * buffer_size + j] + file_size - 1) % file_size];
            fwrite(char_buffer, 1, buffer_size, destination);
        }
        else {
            for (int j = 0; j < file_size % buffer_size; ++j)
                char_buffer[j] = original[(suffix_array[i * buffer_size + j] + file_size - 1) % file_size];
            fwrite(char_buffer, 1, file_size % buffer_size, destination);
        }
    }
    free(char_buffer);
    free(original);
    free(is_small);
    free(is_bucket_end);
    free(suffix_array);
    fclose(destination);
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
