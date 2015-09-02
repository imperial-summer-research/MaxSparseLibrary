
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <queue>

// -- C Library
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <cstdint>

#include "Maxfiles.h"
#include <MaxSLiCInterface.h>

using namespace std;

typedef uint32_t index_t;
typedef float    value_t;
typedef bool     bool_t;
typedef double	 value_rom_t;

class SparseBCSRMatrix {
public:
    // this constructor only apply those design configure variables.
    SparseBCSRMatrix(int r, int c, int depth): r_(r), c_(c), depth_(depth) {
        printf("Design:\n\t%d X %d with %d ROM\n", r_, c_, depth_);
    }

    void MakeDense(int m, int n);
    void ConvertToBCSR();
private:
    // design parameters
    int r_, c_, depth_;
    // matrix parameters
    int m_, n_, nnz_;
    // matrix data
    index_t *rows_, *cols_;
    value_t *vals_;
};

int main(int argc, char *argv[]) {
    int iter = 100;
    int m, n, nnz, r, c, freq;
    // it's a dense matrix
    r = SparseBCSR_R;
    c = SparseBCSR_C;
    n = SparseBCSR_depth;
    m = n;
    nnz = m * n;
    freq = SparseBCSR_freq;

    int loop_length = SparseBCSR_loopLength;

    //SparseBCSRMatrix BCSR_matrix(SparseBCSR_R, SparseBCSR_C, SparseBCSR_depth);
    //BCSR_matrix.MakeDense(m, n);
    //BCSR_matrix.ConvertToBCSR();
//
    printf("%4d X %4d block for %4d X %4d matrix\n", r, c, m, n);
    printf("Number of Non zeros:\t%d\n", nnz);
    printf("Memory size for matrix:\t%.f GB\n", (double)sizeof(index_t) * nnz * 3 / (1024 * 1024 * 1024));

    // Original COO format
    index_t *row = (index_t *) malloc(sizeof(index_t) * nnz);
    index_t *col = (index_t *) malloc(sizeof(index_t) * nnz);
    value_t *val = (value_t *) malloc(sizeof(value_t) * nnz);

    for (int i = 0; i < nnz; i++) {
        row[i] = i / n;
        col[i] = i % n;
        val[i] = 1;
    }
    printf("Generated COO Matrix\n");
    
    // BCSR format: 
    // param{index}: matrix element's column index
    // param{value}: matrix element's value
    // param{input}: vector value input
    vector< vector<index_t> > index_queues(m, vector<index_t>(0));
    vector< vector<value_t> > value_queues(m, vector<value_t>(0));
    for (int i = 0; i < nnz; i++) {
        index_queues[row[i]].push_back(col[i]);
        value_queues[row[i]].push_back(val[i]);
    }
    printf("Transformed to intermediate queues\n");
    vector<index_t> index;
    vector<value_t> value;
    vector<index_t> start;
    // padding
    long long pad_num_row = ceil((double)m/r) * r;
    long long pad_num_col = ceil((double)n/c) * c;
    long long pad_nnz = 0;

    value_t *vec = (value_t *) malloc(sizeof(value_t) * pad_num_col);
    value_t *res = (value_t *) malloc(sizeof(value_t) * pad_num_row);
    for (int i = 0; i < pad_num_col; i++) 
        vec[i] = 1;
    for (int i = 0; i < pad_num_row; i++)
        res[i] = 0.0;

    // generalized iteration
    for (int i = 0; i < pad_num_row; i += r) {
        for (int j = 0; j < pad_num_col; j += c) {
            start.push_back((j == 0));
            bool is_empty = true;
            for (int _i = 0; _i < r; _i++) {
                int _r = i + _i;
                int _c = j;
                if (_r < m && _c < index_queues[_r].size())
                    is_empty = false; 
            }
            if (is_empty) 
                break;
            // for each block
            for (int _j = 0; _j < c; _j ++) {
                for (int _i = 0; _i < r; _i ++) {
                    int _r = i + _i;
                    int _c = j + _j;
                    index_t _index = (_r >= m || _c >= index_queues[_r].size()) ? 0   : index_queues[_r][_c];
                    value_t _value = (_r >= m || _c >= value_queues[_r].size()) ? 0.0 : value_queues[_r][_c];
                    index.push_back(_index);
                    value.push_back(_value);
                    pad_nnz ++;
                }
            }
        }
    }

    // reorder to CSlow sequence
    value_t *value_stream = (value_t *) malloc(sizeof(value_t) * pad_nnz);
    index_t *index_stream = (index_t *) malloc(sizeof(index_t) * pad_nnz);
    index_t *start_stream = (index_t *) malloc(sizeof(index_t) * pad_nnz / (r * c));
    value_t *output_stream = (value_t *) malloc(sizeof(value_t) * pad_nnz / c);

    // several parameters
    int blk_size = r * c;
    int num_blk = pad_nnz / blk_size;
    int arr_size = r * pad_num_col;
    int num_arr = pad_nnz / arr_size;
    int blk_per_arr = num_blk / num_arr;
    int grp_size = r * pad_num_col * loop_length;
    int num_grp = pad_nnz / grp_size;
    int blk_per_grp = num_blk / num_grp;
    int arr_per_grp = num_arr / num_grp;

    printf("Loop length:\t\t %d\n", loop_length);
    printf("Array Size:\t\t %d\n", arr_size);
    printf("Group Size:\t\t %d\n", grp_size);
    printf("Block per Array:\t %d\n", blk_per_arr);
    printf("Block per Group:\t %d\n", blk_per_grp);
    for (int i = 0; i < num_grp; i++) {
        for (int j = 0; j < blk_per_arr; j++) {
            for (int k = 0; k < loop_length; k++) {
                int old_blk_id = i * blk_per_grp + k * blk_per_arr + j;
                int new_blk_id = i * blk_per_grp + j * loop_length + k;
                int new_start_id = new_blk_id * blk_size;
                int old_start_id = old_blk_id * blk_size;
                
                start_stream[new_blk_id] = start[old_blk_id];
                for (int y = 0; y < c; y++)
                    for (int x = 0; x < r; x++) {
                        int new_id = new_start_id + y * r + x;
                        int old_id = old_start_id + y * r + x;
                        value_stream[new_id] = value[old_id];
                        index_stream[new_id] = index[old_id];
                    }
            }
        }
    }

    //for (int i = 0; i < num_blk; i++)
    //    printf("start[%3d] = %d\n", i, start_stream[i]);

    float num_fp_ops = 0;
    num_fp_ops += 1; // multiplier for each nnz
    num_fp_ops += (1-pow(0.5, log2((double)c)+1))*2 - 1 + (double)1/c; // adders for each element
    // SEND TRANSFORMED RESULT TO BOARD
    cout << "Ready to send data to Max3 Board ..."  << endl;
    cout << "Number of non zeros(padded):\t"        << pad_nnz << endl;
    cout << "Number of rows(padded):\t\t"           << pad_num_row << endl;
    cout << "Number of cols(padded):\t\t"           << pad_num_col << endl;
    cout << "Number of blocks(r X c):\t"            << pad_nnz / (r * c) << endl;
    cout << "Frequency for one tick: \t"            << freq << " MHZ" << endl; 
    cout << "Number of fp operations: \t"           << num_fp_ops << endl;

    value_rom_t *rom = (value_rom_t*) malloc(sizeof(value_rom_t) * SparseBCSR_depth);
    for (int i = 0; i < SparseBCSR_depth; i++)
        rom[i] = (i >= pad_num_col) ? 0.0 : (value_rom_t) vec[i];

    //copy(value.begin(), value.end(), value_stream);
    //copy(index.begin(), index.end(), index_stream);
    //copy(start.begin(), start.end(), start_stream);

    cout << "Heating ..." << endl;
    SparseBCSR(pad_nnz, index_stream, start_stream, value_stream, output_stream, rom);

    cout << "Running DFE ..." << endl;
    struct timeval t0, t1;
    gettimeofday(&t0, 0);
    for (int i = 0; i < iter; i++)
        SparseBCSR(pad_nnz, index_stream, start_stream, value_stream, output_stream, rom);
    gettimeofday(&t1, 0);

    cout << "Finished" << endl; 
    double duration = (double)(t1.tv_sec-t0.tv_sec)+(double)(t1.tv_usec-t0.tv_usec)/1e6;
    duration /= iter;
    printf("Total time %lf s per element time %lf us GFLOPS %.6f Eff. GFLOPS %.6f Freq %.2f MHz Bandwidth %.2f MB/s\n", 
        duration, 
        duration * 1e6 / nnz,
        2.0 * nnz / duration / 1e9,
        num_fp_ops * nnz / duration / 1e9,
        1 / (duration / num_blk) / 1e6,
        sizeof(value_t) * nnz / duration / 1e6);

    for (int i = 0; i < num_grp; i++) {
        for (int j = 0; j < loop_length; j++) {
            for (int k = 0; k < r; k++) {
                int row = i * r * loop_length + j * r + k;
                int end_blk_id = i * blk_per_grp + blk_per_grp - loop_length + j;
                int idx = end_blk_id * r + k;
                //printf("row: %d idx: %d\n", row, idx);
                res[row] = output_stream[idx];
            }
        }
    }
    
    value_t *expected = (value_t *) malloc(sizeof(value_t) * m);
    memset(expected, 0, sizeof(value_t) * m);
    for (int i = 0; i < nnz; i++)
        expected[row[i]] += val[i] * vec[col[i]];

    for (int i = 0; i < m; i++)
        if (abs(res[i]-expected[i])/expected[i] > 1e-4) {
            printf("ERROR: [%6d] %15.6f %15.6f\n", i, res[i], expected[i]);
            exit(1);
        }

    printf("OK!\n");
    
    return 0;
}

void SparseBCSRMatrix::MakeDense(int m, int n) {
    //m_ = m, n_ = n, nnz_ = m * n;
//
    //printf("Dense Matrix:\n\t%d X %d = %d\n", m_, n_, nnz_);
    //rows_ = new index_t [nnz_];
    //cols_ = new index_t [nnz_];
    //vals_ = new value_t [nnz_];
//
    //for (int i = 0; i < nnz_; i++) {
    //    rows_[i] = i / n_;
    //    cols_[i] = i % n_;
    //    vals_[i] = (value_t) rand() / RAND_MAX;
    //}
}

void SparseBCSRMatrix::ConvertToBCSR() {
    // BCSR format: 
    // param{index}: matrix element's column index
    // param{value}: matrix element's value
    // param{input}: vector value input
    //vector< vector<index_t> > index_queues(m_, vector<index_t>(0));
    //vector< vector<value_t> > value_queues(m_, vector<value_t>(0));
    //for (int i = 0; i < nnz_; i++) {
    //    index_queues[rows_[i]].push_back(col_[i]);
    //    value_queues[rows_[i]].push_back(val_[i]);
    //}
    //printf("Transformed to intermediate queues\n");
    //vector<index_t> index;
    //vector<value_t> value;
    //vector<index_t> start;
    //// padding
    //long long pad_num_row = ceil((double)m/r) * r;
    //long long pad_num_col = ceil((double)n/c) * c;
    //long long pad_nnz = 0;
//
    //value_t *vec = (value_t *) malloc(sizeof(value_t) * pad_num_col);
    //value_t *res = (value_t *) malloc(sizeof(value_t) * pad_num_row);
    //for (int i = 0; i < pad_num_col; i++) 
    //    vec[i] = 0.01;
    //for (int i = 0; i < pad_num_row; i++)
    //    res[i] = 0.0;
//
    //// generalized iteration
    //for (int i = 0; i < pad_num_row; i += r) {
    //    for (int j = 0; j < pad_num_col; j += c) {
    //        start.push_back((j == 0));
    //        bool is_empty = true;
    //        for (int _i = 0; _i < r; _i++) {
    //            int _r = i + _i;
    //            int _c = j;
    //            if (_r < m && _c < index_queues[_r].size())
    //                is_empty = false; 
    //        }
    //        if (is_empty) 
    //            break;
    //        // for each block
    //        for (int _j = 0; _j < c; _j ++) {
    //            for (int _i = 0; _i < r; _i ++) {
    //                int _r = i + _i;
    //                int _c = j + _j;
    //                index_t _index = (_r >= m || _c >= index_queues[_r].size()) ? 0   : index_queues[_r][_c];
    //                value_t _value = (_r >= m || _c >= value_queues[_r].size()) ? 0.0 : value_queues[_r][_c];
    //                index.push_back(_index);
    //                value.push_back(_value);
    //                pad_nnz ++;
    //            }
    //        }
    //    }
    //}
}