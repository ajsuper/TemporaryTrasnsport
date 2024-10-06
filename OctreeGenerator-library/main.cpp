#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <stdint.h>
#include <bitset>
#include <math.h>
#include <cstdint>
#include <algorithm>
#include <ctime>
#include <omp.h>

std::vector<std::vector<std::vector<int>>> createVoxelGrid(int size){
    std::vector<std::vector<std::vector<int>>> voxels(size, std::vector<std::vector<int>>(size, std::vector<int>(size, 0)));
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            for(int k = 0; k < size; k++){
                voxels[i][j][k] = 0;
            }
        }
    }
    return voxels;
}
//write a vector of uint8_t's to a file
void writeUint32Vector(std::vector<uint32_t> vect, std::string fileDirectory){
    std::ofstream outFile(fileDirectory, std::ios::binary);
    if(!outFile){
        std::cout << "ERROR in 'writeUint8Vector': Failed to open " << fileDirectory << std::endl;
        return;
    }
    size_t size = vect.size();
    outFile.write(reinterpret_cast<const char*>(&size), sizeof(size));
    outFile.write(reinterpret_cast<const char*>(vect.data()), vect.size() * sizeof(uint32_t));
    outFile.close();
}

//read a file containing a vector of uint8_t's
std::vector<uint32_t> readUint32Vector(std::string fileDirectory){
    std::ifstream inFile(fileDirectory, std::ios::binary);  
    if (!inFile) {
        std::cerr << "Error opening file for reading." << std::endl;
    }
    size_t size;
    inFile.read(reinterpret_cast<char*>(&size), sizeof(size));

    std::vector<uint32_t> readNumbers(size);

    inFile.read(reinterpret_cast<char*>(readNumbers.data()), size * sizeof(uint32_t));
    inFile.close();

    //std::cout << "Read data: " << std::endl;
    for (uint32_t num : readNumbers) {
        //std::cout << std::bitset<32>(num) << std::endl;
    } 
    return readNumbers; 
}

int voxInBounds(std::vector<std::vector<std::vector<int>>>& voxels, int wholeVoxelResolution, const int octantPosition[3], const int octantPositionExtents[3]){
    int voxCountInBounds = 0;
    for (int x = 0; x < 3; x++){
        if(octantPositionExtents[x] > wholeVoxelResolution || octantPosition[x] < 0){
            return 0;
        }
    }
    for(int x = octantPosition[0]; x < octantPositionExtents[0]; x++){
        for(int y = octantPosition[1]; y < octantPositionExtents[1]; y++){
            for(int z = octantPosition[2]; z < octantPositionExtents[2]; z++){
                if(voxels[x][y][z] != 0){
                    voxCountInBounds += 1;
                }
                //voxels[x][y][z] = 1;
            }
        }
    }

    return voxCountInBounds;
}

uint16_t buildMask(std::vector<std::vector<std::vector<int>>>& voxels, const int& checkX, const int& checkY, const int& checkZ, const int& depthResolution, int wholeVoxelResolution){
    uint16_t masks = 0;
    int childSize = int(depthResolution/2);
    int x = 0;
    int y = 1;
    int z = 2;

    //std::cout << checkX << " " << checkY << " " << checkZ << " " << depthResolution << std::endl;

    int octants[8][3] = {
        {checkX, checkY, checkZ},
        {checkX, checkY, checkZ+childSize},
        {checkX, checkY+childSize, checkZ},
        {checkX, checkY+childSize, checkZ+childSize},
        {checkX+childSize, checkY, checkZ},
        {checkX+childSize, checkY, checkZ+childSize},
        {checkX+childSize, checkY+childSize, checkZ},
        {checkX+childSize, checkY+childSize, checkZ+childSize}
    };

    

    for(int i = 0; i < 8; i++){
        //std::cout << octants[i][0] << " " << octants[i][1] << " " << octants[i][2] << " " << childSize << std::endl;
    }
    if(childSize > 1){
        for(int i = 0; i < 8; i++){
            int octantPosition[3] = {octants[i][x], octants[i][y], octants[i][z]};
            int octantPositionExtents[3] = {octants[i][x]+childSize, octants[i][y]+childSize, octants[i][z]+childSize};
            
            if(voxInBounds(voxels, wholeVoxelResolution, octantPosition, octantPositionExtents) > 0){
                //std::cout << octantPosition[0] << " " << octantPosition[1] << " " << octantPosition[2] << std::endl;
                masks |= (1 << (15-i));
            }
        }
    }

    if(childSize == 1){
        for(int i = 0; i < 8; i++){
            int octantPosition[3] = {octants[i][x], octants[i][y], octants[i][z]};
            int octantPositionExtents[3] = {octants[i][x]+childSize, octants[i][y]+childSize, octants[i][z]+childSize};
            if(voxInBounds(voxels, wholeVoxelResolution, octantPosition, octantPositionExtents) > 0){
                masks = masks | (1 << (15-i));
                masks = masks | (1 << (7-i));
            }
        }
    }

    return masks;
}

uint64_t createZOrderIndex(const int& x, const int& y, const int& z, const int& depth) {
    uint64_t z_order = 0;

    for (int i = 0; i < depth; ++i) {
        // Swapping x and y axes
        uint64_t bitY = (z >> i) & 1;  // Extract the i-th bit of y
        uint64_t bitX = (y >> i) & 1;  // Extract the i-th bit of x
        uint64_t bitZ = (x >> i) & 1;  // Extract the i-th bit of z

        // Interleave them with swapped x and y (y -> x -> z)
        z_order |= (bitY << (3 * i)) | (bitX << (3 * i + 1)) | (bitZ << (3 * i + 2));
    }

    return z_order;
}


void sortCoordinatesToZOrder(std::vector<std::array<int, 3>>& coordinates, uint32_t depth) {
    std::sort(coordinates.begin(), coordinates.end(), [depth](const auto& a, const auto& b) {
        // Use Z-order index for sorting the coordinates.
        return createZOrderIndex(a[0], a[1], a[2], depth) < 
               createZOrderIndex(b[0], b[1], b[2], depth);
    });
}

//Build all of the masks for each node at the current depth.
//Defined by the depth in the tree and the voxel width count.
std::vector<uint16_t> buildMasksForWholeDepth(std::vector<std::vector<std::vector<int>>>& voxels, int depthInOctree, int voxelWholeResolution){
    int voxelDepthResolution = pow(2, depthInOctree);
    int wholeResToDepthResRatio = voxelWholeResolution/voxelDepthResolution;
    int index = 0;

    std::vector<uint16_t> thisDepthMasks;
    std::vector<std::array<int, 3>> coordinates;
    coordinates.resize(voxelDepthResolution*voxelDepthResolution*voxelDepthResolution);
    
    for (int i = 0; i < voxelDepthResolution; i++) {
        for (int j = 0; j < voxelDepthResolution; j++) {
            for (int k = 0; k < voxelDepthResolution; k++){
                coordinates[index] = {i, j, k};
                index += 1;
            }
        }
    }
    int bitLength = log2(voxelWholeResolution); // Number of bits for the coordinate
    sortCoordinatesToZOrder(coordinates, bitLength);
    //Loop over each of the coordinates and create a mask then add that to the whole list of masks.
    for(int i = 0; i < coordinates.size(); i++){
        std::array<int, 3> coord = coordinates[i];
        
        uint16_t temporaryMask = buildMask(voxels, wholeResToDepthResRatio*coord[0], wholeResToDepthResRatio*coord[1], wholeResToDepthResRatio*coord[2], wholeResToDepthResRatio, voxelWholeResolution);
        //std::cout << std::bitset<16>(temporaryMask) << " " << coord[0]*wholeResToDepthResRatio << " " << coord[1]*wholeResToDepthResRatio << " " << coord[2]*wholeResToDepthResRatio << " " << std::endl;
        if(temporaryMask != 0){
            thisDepthMasks.push_back(temporaryMask);
        }
    }
    return thisDepthMasks;
}



std::vector<uint32_t> addPointers(std::vector<uint32_t>& octree){
    uint16_t vldMaskMask = 0b1111111100000000;
    uint16_t lefMaskMask = 0b0000000011111111;
    int parentCount = 0;
    int childCounter = 0;
    for(int address = 0; address < octree.size(); address++){
        uint32_t* current = &octree[address];
        uint16_t childPointer;
        uint8_t validMask = (*current & vldMaskMask) >> 8;
        uint8_t leafMask = (*current & lefMaskMask) >> 0;
        if(validMask != 0 && leafMask == 0){ 
            childPointer = childCounter-address+1;
            parentCount += 1;
            for(uint8_t i = 0; i < 8; i++){
                if((validMask & (1 << i)) != 0){
                    childCounter += 1;
                }
            }
        } else {
            childPointer = 0;
        }
        octree[address] = *current | childPointer << 16;
    }
    return octree;
}
//Create a SVO heirarchy based off of a voxel grid of a size that is a power of 2
std::vector<uint32_t> createOctree(std::vector<std::vector<std::vector<int>>>& voxels, int voxelWholeResolution){
    std::vector<uint32_t> octree;
    int levelsOfDepth = int(log2(voxelWholeResolution));
    for(int depthInOctree = 0; depthInOctree < levelsOfDepth; depthInOctree++){
        std::vector<uint16_t> slicesAtDepth = buildMasksForWholeDepth(voxels, depthInOctree, voxelWholeResolution);
        int sliceCountAtDepth = slicesAtDepth.size();
        for(int j = 0; j < sliceCountAtDepth; j++){
            octree.push_back(slicesAtDepth[j]);
        }
    }
    octree = addPointers(octree);
    return octree;
}


int main(){
    readUint32Vector("test8x8x8.bin");
    int size = 256;
    int voxelCount = 0;
    std::vector<std::vector<std::vector<int>>> voxels = createVoxelGrid(size);
    for(int x = 0; x < size; x++){
        for(int y = 0; y < size; y++){
            for(int z = 0; z < size; z++){ 

                int dx = x - size;
                int dy = y - size/2;
                int dz = z - size/2;
                float distanceSquared = dx * dx + dy * dy + dz * dz;
                if(distanceSquared <= (size/2)*(size/2) && distanceSquared >= ((size-2)/2)*((size-2)/2)){
                //if(distanceSquared <= (size/2)*(size/2)){
                    voxels[x][y][z] = 1;
                    voxelCount += 1;
                }
                
           
                if(y < x/2){
                    //voxels[x][y][z] = 1;
                }
                if(y == 2){
                    //voxels[x][y][z] = 1;
                }
                if(y < sin(float(x)/5)*13+25){
                    //voxels[x][y][z] = 1;
                }
                //voxels[x][y][z] = 1;
            }
        }
           
    }

    //
    //voxels[15][15][15] = 1;
    //voxels[0][4][0] = 1;
    //voxels[0][0][1] = 1;
    //voxels[1][0][2] = 1;

    std::vector<uint32_t> octree;
    std::clock_t start = std::clock();
    octree = createOctree(voxels, size);
    std::clock_t end = std::clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "time: " << elapsed*1000 << " ms\n" << std::endl;
    std::cout << "voxelCount:" << voxelCount << std::endl;
    uint64_t largest = 0;
    for(int i = 0; i < octree.size(); i++){
        if(octree[i] > largest){
            largest = octree[i];
        }
        //std::cout << std::bitset<16>(octree[i] >> 16) << "|" << std::bitset<16>(octree[i]) << std::endl;
        //std::cout << "p:" << std::bitset<16>(octree[i] >> 16) << "m:" << std::bitset<16>(octree[i]) << std::endl;
    }
    std::cout << "largest" << std::bitset<32>(largest) << std::endl;

    uint16_t masktest2 = buildMask(voxels, 0, 0, 0, 4, voxels.size());
    int bound1[3] = {4,0,8};
    int bound2[3] = {8,4,12};
    bool jeff = voxInBounds(voxels, voxels.size(), bound1, bound2);
    std::cout << std::bitset<16>(masktest2) << std::endl;
    writeUint32Vector(octree, "test8x8x8.bin");

    
    return 0;
}

