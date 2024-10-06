#include </usr/include/c++/13/fstream>
#include </usr/include/c++/13/sstream>
#include </usr/include/c++/13/iostream>
#include </usr/include/GL/glew.h>
#include </usr/include/GL/freeglut.h>
#include </usr/include/c++/13/math.h>
#include </usr/include/c++/13/vector>
#include <bitset>
#include <chrono>


auto last = std::chrono::high_resolution_clock::now();
auto now = std::chrono::high_resolution_clock::now();
//Initialize our shaderProgram
GLuint shaderProgram;

float timeValue = 0.0f;
float movementSpeed = 0.3f;
struct camera {
    float position[3] = {33, 10, 36};
    float direction[3] = {0, 0, -1};
    float up[3] = { 0, 1, 0 };
} cam;

struct boxStruct {
    float position[3];
    float size;
};


struct recursiveBoxStruct {
    int address;
    float position[3];
    float size;
    int steps;
};


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
        std::cout << std::bitset<32>(num) << std::endl;
    } 
    return readNumbers; 
}


std::vector<boxStruct> voxelsToBoxes(std::vector<std::vector<std::vector<int>>> voxels){
    std::vector<boxStruct> boxes;
    for(int x = 0; x < voxels.size(); x++){
        for(int y = 0; y < voxels.size(); y++){
            for(int z = 0; z < voxels.size(); z++){
                if(voxels[x][y][z] != 0){
                    boxStruct box;
                    box.position[0] = x;
                    box.position[1] = y;
                    box.position[2] = z;
                    box.size = 1.0f;
                    boxes.push_back(box);
                }
            }
        }
    }
    return boxes;
}


void setVoxels(std::vector<std::vector<std::vector<int>>>& voxels){
    for(int x = 0; x < voxels.size(); x++){
        for(int y = 0; y < voxels.size(); y++){
            for(int z = 0; z < voxels.size(); z++){
                int dx = x - 2;
                int dy = y - 2;
                int dz = z - 2;
                float distanceSquared = dx * dx + dy * dy + dz * dz;
                if(distanceSquared <= 2*2){
                    voxels[x][y][z] = 1;
                }
            }
        }
    }
    return;
}


void passBoxesToFrag(std::vector<boxStruct> boxes){
    GLuint ssbo;
    glGenBuffers(1, &ssbo);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
    glBufferData(GL_SHADER_STORAGE_BUFFER, boxes.size() * sizeof(boxStruct), boxes.data(), GL_STATIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssbo);
}

void passOctreeToFrag(std::vector<uint32_t> octree){
    GLuint ssbo;
    glGenBuffers(1, &ssbo);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
    glBufferData(GL_SHADER_STORAGE_BUFFER, octree.size() * sizeof(uint32_t), octree.data(), GL_STATIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssbo);
}

std::vector<boxStruct> octreeToBoxes(std::vector<uint32_t> octree, int levelsOfDepth){
    std::vector<recursiveBoxStruct> stack;
    std::vector<boxStruct> boxList;
    
    recursiveBoxStruct rootBox;
    boxStruct rootBoxReal;
    
    // Initialize the root box and push onto the stack
    rootBox.address = 0;
    rootBox.position[0] = 5;
    rootBox.position[1] = 0;
    rootBox.position[2] = 0;
    rootBox.size = 2.0;
    rootBox.steps = 0;
    
    rootBoxReal.position[0] = rootBox.position[0];
    rootBoxReal.position[1] = rootBox.position[1];
    rootBoxReal.position[2] = rootBox.position[2];
    rootBoxReal.size = rootBox.size;
    
    boxList.push_back(rootBoxReal);
    stack.push_back(rootBox);
    
    while(stack.size() > 0 && octree.size() > 0){
        recursiveBoxStruct recursiveBox = stack.back();
        stack.pop_back();
        
        // Stop if depth exceeds the level of the octree
        if(recursiveBox.steps > levelsOfDepth){
            return boxList;
        }
        
        uint32_t node = octree[recursiveBox.address];
        uint8_t validMask = (node >> 8) & 0xFF;
        uint8_t leafMask = node & 0xFF;
        uint16_t pointer = (node >> 16);
        int firstChildAddress = pointer + recursiveBox.address;
        int childCount = 0;

        for(int octant = 0; octant < 8; octant++){
            int octantArray[3] = {0, 0, 0};
            if(octant == 1){ octantArray[0] = 1; octantArray[1] = 0; octantArray[2] = 0; };
            if(octant == 2){ octantArray[0] = 0; octantArray[1] = 1; octantArray[2] = 0; };
            if(octant == 3){ octantArray[0] = 1; octantArray[1] = 1; octantArray[2] = 0; };
            if(octant == 4){ octantArray[0] = 0; octantArray[1] = 0; octantArray[2] = 1; };
            if(octant == 5){ octantArray[0] = 1; octantArray[1] = 0; octantArray[2] = 1; };
            if(octant == 6){ octantArray[0] = 0; octantArray[1] = 1; octantArray[2] = 1; };
            if(octant == 7){ octantArray[0] = 1; octantArray[1] = 1; octantArray[2] = 1; };
            
            float position[3] = {
                recursiveBox.position[0] + (octantArray[0] * recursiveBox.size / 2),
                recursiveBox.position[1] + (octantArray[1] * recursiveBox.size / 2),
                recursiveBox.position[2] + (octantArray[2] * recursiveBox.size / 2)
            };

            // Fix operator precedence
            if((validMask & (1 << (7 - octant))) != 0){
                int thisAddress = firstChildAddress + childCount;
                childCount += 1;
                
                boxStruct newBox;
                newBox.position[0] = position[0];
                newBox.position[1] = position[1];
                newBox.position[2] = position[2];
                newBox.size = recursiveBox.size / 2;
                
                
                boxList.push_back(newBox);
                
                
                
                

                
                // Fix operator precedence and add recursive box if it's not a leaf
                if((leafMask & (1 << (7 - octant))) == 0){
                    recursiveBoxStruct newBoxRecursive;
                    newBoxRecursive.address = thisAddress;
                    newBoxRecursive.position[0] = position[0];
                    newBoxRecursive.position[1] = position[1];
                    newBoxRecursive.position[2] = position[2];
                    newBoxRecursive.size = recursiveBox.size / 2;
                    newBoxRecursive.steps = recursiveBox.steps + 1;  // Increment steps
                    
                    stack.push_back(newBoxRecursive);
                }
            }
        }
    }
    
    return boxList;
}



void defineTriangles() {
    glBegin(GL_TRIANGLES);
    glVertex3f(-1.0f, -1.0f, 0.0f);
    glVertex3f(1.0f, -1.0f, 0.0f);
    glVertex3f(-1.0f, 1.0f, 0.0f);
    glVertex3f(1.0f, -1.0f, 0.0f);
    glVertex3f(-1.0f, 1.0f, 0.0f);
    glVertex3f(1.0f, 1.0f, 0.0f);
    glEnd();
}


void reinitializeEachFrame() {
    glUseProgram(shaderProgram);
    timeValue = timeValue + 0.05;
}


void passValues() {
    GLint resolutionLoc = glGetUniformLocation(shaderProgram, "resolution");
    float res[2] = { static_cast<GLfloat>(glutGet(GLUT_WINDOW_WIDTH)), static_cast<GLfloat>(glutGet(GLUT_WINDOW_HEIGHT)) };
    glUniform2fv(resolutionLoc, 1, res);

    GLint timeLocation = glGetUniformLocation(shaderProgram, "time");
    if (timeLocation != -1) {
        glUniform1f(timeLocation, (timeValue / 2000));
    }

    GLint cameraPosLoc = glGetUniformLocation(shaderProgram, "cameraPos");
    GLint cameraDirLoc = glGetUniformLocation(shaderProgram, "cameraDir");
    glUniform3fv(cameraPosLoc, 1, cam.position);
    glUniform3fv(cameraDirLoc, 1, cam.direction);
    
    //Create the octree texture

    //glTexImage1D(GL_TEXTURE_1D, 0, GL_R32F, textureWidth, 0, GL_RED, GL_FLOAT, textureData);

}

 
void keys(unsigned char key, int x, int y) {
    float right[3];
    // Cross product
    right[0] = cam.up[1] * cam.direction[2] - cam.up[2] * cam.direction[1];
    right[1] = cam.up[2] * cam.direction[0] - cam.up[0] * cam.direction[2];
    right[2] = cam.up[0] * cam.direction[1] - cam.up[1] * cam.direction[0];
    // Normalize
    float norm = sqrt(right[0] * right[0] + right[1] * right[1] + right[2] * right[2]);
    right[0] /= norm;
    right[1] /= norm;
    right[2] /= norm;
    float up[3];
    up[1] = 1;
    switch (key) {
        //Direction
        case 'e':
            cam.direction[0] = cam.direction[0] * cos(-movementSpeed) - cam.direction[2] * sin(-movementSpeed);
            cam.direction[2] = cam.direction[0] * sin(-movementSpeed) + cam.direction[2] * cos(-movementSpeed);
            break;
        case 'q':
            cam.direction[0] = cam.direction[0] * cos(movementSpeed) - cam.direction[2] * sin(movementSpeed);
            cam.direction[2] = cam.direction[0] * sin(movementSpeed) + cam.direction[2] * cos(movementSpeed);
            break;
        // Position
        case 'w':
            cam.position[0] += cam.direction[0] * movementSpeed;
            cam.position[1] += cam.direction[1] * movementSpeed;
            cam.position[2] += cam.direction[2] * movementSpeed;
            break;
        case 'd':
            cam.position[0] += right[0] * movementSpeed;
            cam.position[1] += right[1] * movementSpeed;
            cam.position[2] += right[2] * movementSpeed;
            break;
        case 's':
            cam.position[0] -= cam.direction[0] * movementSpeed;
            cam.position[1] -= cam.direction[1] * movementSpeed;
            cam.position[2] -= cam.direction[2] * movementSpeed;
            break;
        case 'a':
            cam.position[0] -= right[0] * movementSpeed;
            cam.position[1] -= right[1] * movementSpeed;
            cam.position[2] -= right[2] * movementSpeed;
            break;    
        case 'r':
            cam.position[1] += up[1]*movementSpeed;  
            break;
        case 'f':
            cam.position[1] -= up[1]*movementSpeed;
    }
}


std::string loadShaderSource(const char* filepath) {
    std::ifstream shaderFile;
    shaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
        shaderFile.open(filepath);
        std::stringstream shaderStream;
        shaderStream << shaderFile.rdbuf();
        shaderFile.close();
        return shaderStream.str();
    }
    catch (std::ifstream::failure& e) {
        std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ" << std::endl;
        return "";
    }
}


void compileAndUseProgramShader() {
    //Get shader code
    std::string vertexCode = loadShaderSource(R"(./shaders/vertexShader.glsl)");
    std::string fragmentCode = loadShaderSource(R"(./shaders/fragmentShader.frag)");

    //Convert it to a string
    const char* vertexShaderSource = vertexCode.c_str();
    const char* fragmentShaderSource = fragmentCode.c_str();

    // Define the vertexShaderSource as a GL Vertex Shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    // Define the fragmentShaderSource as a GL Fragment Shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    // Combine the seperate shaders into the one shader shaderProgram
    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    // Now that we have shaderProgram we no longer need each seperate shader so we can delete them
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
}


void display(void)
{
    static float last_time = 0;
    static float fps_start = 0;

    now = std::chrono::high_resolution_clock::now();
    auto totalTimeForFrame = now-last; 
    auto durationInMilliseconds = std::chrono::duration_cast<std::chrono::microseconds>(totalTimeForFrame);
    std::cout << durationInMilliseconds.count()/ 1000.0 << " ms" << std::endl;
    last = now;

    reinitializeEachFrame();
    passValues();
    defineTriangles();
    glFinish();
    glutPostRedisplay();
}


int main(int argc, char* argv[]) {
    std::vector<uint32_t> octree = readUint32Vector("/home/andrew/vscode/cpp/OctreeGenerator-library/test8x8x8.bin");
    const int size = 64;
    std::vector<std::vector<std::vector<int>>> voxels(size, std::vector<std::vector<int>>(size, std::vector<int>(size)));
    std::vector<boxStruct> boxesFromOctree = octreeToBoxes(octree, log2(size));
    glutInit(&argc, argv);
    glutInitWindowPosition(1, 1);
    glutInitWindowSize(1900, 1900);
    glutCreateWindow("RayTracer!");
    glewInit();
    compileAndUseProgramShader();

    //setVoxels(voxels);
    voxels[0][0][0] = 1;

    passBoxesToFrag(boxesFromOctree);
    passOctreeToFrag(octree);
    
    glutDisplayFunc(display);
    glutKeyboardFunc(keys);
    glutMainLoop();
    return 0;
}