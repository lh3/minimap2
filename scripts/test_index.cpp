#include <iostream>


int main(int argc, const char** argv) {
    
    unsigned blockdim = atoi(argv[1]);
    unsigned range = atoi(argv[2]);
    unsigned i_range = atoi(argv[3]);
    unsigned buffer_size = atoi(argv[4]);

    unsigned start_idx = 3;

    int anchor_offset[blockdim] = {0};
    for (unsigned tid = 0; tid < blockdim; ++tid) {
        anchor_offset[tid] = tid + (start_idx - tid + blockdim) / blockdim * blockdim;
        std::cout << "tid: " << tid << " " << anchor_offset[tid] << std::endl;
    }

    // unsigned anchors[blockdim][buffer_size] = {0};
    // unsigned i = 0;
    // for (unsigned tid = 0; tid < blockdim; ++tid) {
    //     // update when drop becomes 0
    //     // TODO: think about how to update the array of anchors
    //     for (unsigned j = tid+1, index=0; j < i+range+1; j += blockdim, index++) {
    //         anchors[tid][index] = j;
    //     }
    // }
    // // print the anchors
    // for (unsigned tid = 0; tid < blockdim; ++tid) {
    //     std::cout << "tid: " << tid << " [";
    //     for (unsigned j = 0; j < buffer_size; ++j) {
    //         std::cout << anchors[tid][j] << " ";
    //     }
    //     std::cout << "]" << std::endl;
    // }

    for (unsigned i = start_idx; i < i_range; ++i) {
        std::cout << "anchor i: " << i << std::endl;
        for (unsigned tid = 0; tid < blockdim; ++tid) {
            std::cout << "tid: " << tid << ": ";
            for (unsigned j = anchor_offset[tid]; j < i+range+1; j += blockdim) {
                std::cout << j << " ";
            }
            if (anchor_offset[tid] <= i+1) {
                anchor_offset[tid] += blockdim;
            }
            std::cout << std::endl;
            // unsigned drop = (i+blockdim-tid-1) % blockdim;
            // unsigned offset = (i+blockdim-tid-1) / blockdim;
            // // remove old when drop becomes 0
            // // add new when drop becomes (blockdim - range%blockdim)
            // // TODO: think about how to update the array of anchors efficiently
            // // easiest, generate a mask with ballot and use cooperative groups
            // std::cout << "tid: " << tid  << " offset: " << offset << " drop: " << drop << " (";
            // if (!drop) {
            //     anchors[tid][(offset+buffer_size-1)%buffer_size] = 0;
            // }
            // for (unsigned j = offset*blockdim+tid+1, index=0; j < i+range+1; j += blockdim, index++) {
            //     if (drop == blockdim - range%blockdim && index == buffer_size-1) {
            //         anchors[tid][(index+offset)%buffer_size] = j;
            //     }
            //     std::cout << j << " : " << anchors[tid][(index+offset)%buffer_size] << ", ";
            //     // anchors[tid][index] = j;
            // }
            // std::cout << ") " << std::endl;
        }
        std::cout << std::endl;
        // for (unsigned tid = 0; tid < blockdim; ++tid) {
        //     std::cout << "tid: " << tid << " [";
        //     for (unsigned j = 0; j < buffer_size; ++j) {
        //         std::cout << anchors[tid][j] << " ";
        //     }
        //     std::cout << "]" << std::endl;
        // }

    }
    
    
    return 0;
}
