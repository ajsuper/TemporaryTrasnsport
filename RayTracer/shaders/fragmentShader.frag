#version 450 core
out vec4 FragColor;
in vec2 TexCoords;
uniform float time;
uniform vec2 resolution;
uniform vec3 cameraPos;
uniform vec3 cameraDir;
const int MAX_STACK_SIZE = 10; // Maximum size for stacks
const int MAX_RAY_STEPS = 30;
vec3 invRayDir;

struct rayStruct {
    vec3 origin;
    vec3 direction;
    vec3 color;
};


struct boxStruct {
    vec3 position;
    float size;
};


layout(std430, binding = 0) buffer BoxesBuffer {
    boxStruct boxes[];
};


layout(std430, binding = 1) buffer OctreeBuffer {
    uint octree[];
};

boxStruct boxStack[MAX_STACK_SIZE];
uint nodeIndexStack[MAX_STACK_SIZE];
uint octant1Stack[MAX_STACK_SIZE];

int boxStackQuantity = 0;
int nodeIndexStackQuantity = 0;
int octant1StackQuantity = 0;


float rand(vec2 uv, float timeX){
    float n = fract(sin(dot(uv*timeX, vec2(12.9898, 78.233))) * 43758.5453);
    return n;
}


vec3 rayStartDirection(vec2 uv, vec2 res, vec3 cameraPosition, vec3 cameraDirection, float fov) {
    float aspectRatio = res.x / res.y;
    float scale = tan(radians(fov * 0.5));
    vec2 pixelNDC = vec2(
        ((-uv.x) * scale) * aspectRatio, // Apply aspect ratio to x-coordinate
        uv.y * scale
    );
    vec3 targetPosition = cameraPosition + cameraDirection;
    vec3 lookat = normalize(targetPosition - cameraPosition);
    vec3 right = cross(lookat, vec3(0, 1, 0));
    vec3 actualUp = cross(right, lookat);
    vec3 rayDirection = normalize(pixelNDC.x * right + pixelNDC.y * actualUp + lookat);
    return rayDirection;
}


vec3 calculatePlaneIntersections(vec3 plane, rayStruct castedRay) {
    if (abs(castedRay.direction.x) < 1e-6) castedRay.direction.x = 1e-6;
    if (abs(castedRay.direction.y) < 1e-6) castedRay.direction.y = 1e-6;
    if (abs(castedRay.direction.z) < 1e-6) castedRay.direction.z = 1e-6;
    
    vec3 intersection;
    intersection.x = (plane.x - castedRay.origin.x) / castedRay.direction.x;
    intersection.y = (plane.y - castedRay.origin.y) / castedRay.direction.y;
    intersection.z = (plane.z - castedRay.origin.z) / castedRay.direction.z;
    return intersection;
}


uint moveOctant1(uint octant1, uint sideToMoveTo) {
    uint newOctant1 = 0;

    if(sideToMoveTo == 1) {
        if(octant1 == 0){ newOctant1 = 4; } else
        if(octant1 == 1){ newOctant1 = 5; } else
        if(octant1 == 2){ newOctant1 = 6; } else
        if(octant1 == 3){ newOctant1 = 7; } else
        if(octant1 == 4){ newOctant1 = 0; } else
        if(octant1 == 5){ newOctant1 = 1; } else
        if(octant1 == 6){ newOctant1 = 2; } else
        if(octant1 == 7){ newOctant1 = 3; }
    } else if(sideToMoveTo == 2) {
        if(octant1 == 0){ newOctant1 = 2; } else
        if(octant1 == 1){ newOctant1 = 3; } else
        if(octant1 == 2){ newOctant1 = 0; } else
        if(octant1 == 3){ newOctant1 = 1; } else
        if(octant1 == 4){ newOctant1 = 6; } else
        if(octant1 == 5){ newOctant1 = 7; } else
        if(octant1 == 6){ newOctant1 = 4; } else
        if(octant1 == 7){ newOctant1 = 5; }
    } else if(sideToMoveTo == 3) {
        if(octant1 == 0){ newOctant1 = 1; } else
        if(octant1 == 1){ newOctant1 = 0; } else
        if(octant1 == 2){ newOctant1 = 3; } else
        if(octant1 == 3){ newOctant1 = 2; } else
        if(octant1 == 4){ newOctant1 = 5; } else
        if(octant1 == 5){ newOctant1 = 4; } else
        if(octant1 == 6){ newOctant1 = 7; } else
        if(octant1 == 7){ newOctant1 = 6; }
    }
    
    return newOctant1;
}


vec2 checkBox(boxStruct box, rayStruct castedRay){
    vec3 invDir = 1.0 / castedRay.direction;
    vec3 t0 = (box.position - castedRay.origin) * invDir;
    vec3 t1 = ((box.position + box.size) - castedRay.origin) * invDir;

    vec3 tsmaller = min(t0, t1);
    vec3 tbigger = max(t0, t1);


    float tmin_candidate = max(max(tsmaller.x, tsmaller.y), tsmaller.z);
    float tmax_candidate = min(min(tbigger.x, tbigger.y), tbigger.z);

    if(tmin_candidate > tmax_candidate || tmax_candidate < 0){
        return vec2(-1, 0);
    }

    uint intersectedSide = 0;
    // Determine which axis contributed to tmin
    if (tmin_candidate == tsmaller.x) {
        intersectedSide = 1; // X-axis
    } else if (tmin_candidate == tsmaller.y) {
        intersectedSide = 2; // Y-axis
    } else {
        intersectedSide = 3; // Z-axis
    }

    float tmin = tmin_candidate;
    float tmax = tmax_candidate;

    
    // Check for intersection
    return vec2(tmin, intersectedSide);
}


vec2 checkBoxInverted(boxStruct box, rayStruct castedRay){
    vec3 invDir = 1.0 / castedRay.direction;
    vec3 t0 = (box.position - castedRay.origin) * invDir;
    vec3 t1 = ((box.position + box.size) - castedRay.origin) * invDir;

    vec3 tsmaller = min(t0, t1);
    vec3 tbigger = max(t0, t1);

    float tmin_candidate = max(max(tsmaller.x, tsmaller.y), tsmaller.z);
    float tmax_candidate = min(min(tbigger.x, tbigger.y), tbigger.z);

    if (tmin_candidate > tmax_candidate || tmax_candidate < 0.0) {
        return vec2(-1.0, 0.0); // No intersection
    }

    uint exitSide = 0u;
    // Determine which axis contributed to tmax
    if (tmax_candidate == tbigger.x) {
        exitSide = 1u; // X-axis
    } else if (tmax_candidate == tbigger.y) {
        exitSide = 2u; // Y-axis
    } else {
        exitSide = 3u; // Z-axis
    }

    // Return tmax and exitSide
    return vec2(tmax_candidate, float(exitSide));
}



uint calculateOctant1BasedOfEdgePos(vec3 rayEdgePos) {
    bool minX = true;
    bool minY = true;
    bool minZ = true;
    
    // Adjust these comparisons to the correct axes
    if(rayEdgePos.x > 0.5) { minX = false; }  // X-axis comparison
    if(rayEdgePos.y > 0.5) { minY = false; }  // Y-axis comparison
    if(rayEdgePos.z > 0.5) { minZ = false; }  // Z-axis comparison
    
    // Now the bits will correspond correctly
    if(minX && minY && minZ) { return 0; } else
    if(minX && minY && !minZ) { return 1; } else
    if(minX && !minY && minZ) { return 2; } else
    if(minX && !minY && !minZ) { return 3; } else
    if(!minX && minY && minZ) { return 4; } else
    if(!minX && minY && !minZ) { return 5; } else
    if(!minX && !minY && minZ) { return 6; } else
    if(!minX && !minY && !minZ) { return 7; } else
    
    return 8; // Error case, shouldn't be reached
}



uint calculateSiblingsBeforeThisOctant1(uint node, uint checkOctant1){ //Verified: Push
    uint siblingsBeforeThisOct = 0;
    uint childCount = 0;
    uint validMask = (node >> 8) & 0x00FF;
    for(uint i = 0; i < 8; i++){
        if(i == checkOctant1){ siblingsBeforeThisOct = childCount; }
        if((validMask & (1 << (7-i))) != 0){ childCount += 1; }
    }
    return siblingsBeforeThisOct;
}


uint rayBoxOctant1(rayStruct ray, boxStruct box){ //Verified: Push
    vec2 intersectDat = checkBox(box, ray);
    vec3 intersectPosition = ray.origin+ray.direction*intersectDat[0];
    vec3 intersectPositionUV = (intersectPosition-box.position)/box.size;
    uint octant1 = calculateOctant1BasedOfEdgePos(intersectPositionUV);
    return octant1;
}

boxStruct newBoxFromParentAndOctant1v2(boxStruct parentBox, uint checkOctant1) {
    float childSize = parentBox.size * 0.5;
    boxStruct childBox;
    vec3 offsetPosition = vec3(0.0, 0.0, 0.0);

    // Ensure that GLSL bitwise operations work correctly
    if ((checkOctant1 & 1u) != 0u) { offsetPosition.x = childSize; }  // X-axis
    if ((checkOctant1 & 2u) != 0u) { offsetPosition.z = childSize; }  // Z-axis
    if ((checkOctant1 & 4u) != 0u) { offsetPosition.y = childSize; }  // Y-axis

    childBox.position = parentBox.position + offsetPosition;
    childBox.size = childSize;
    
    return childBox;
}


boxStruct newBoxFromParentAndOctant1(boxStruct parentBox, uint checkOctant1){ //Verified: Push
    float childSize = parentBox.size*0.5;
    boxStruct childBox;
    vec3 offsetPosition = vec3(0.0, 0.0, 0.0);
    if(checkOctant1 == 0){

    } else if(checkOctant1 == 1){
        offsetPosition.z = childSize;
    } else if(checkOctant1 == 2){
        offsetPosition.y = childSize;
    } else if(checkOctant1 == 3){
        offsetPosition.z = childSize;
        offsetPosition.y = childSize;
    } else if(checkOctant1 == 4){
        offsetPosition.x = childSize;
    } else if(checkOctant1 == 5){
        offsetPosition.z = childSize;
        offsetPosition.x = childSize;
    } else if(checkOctant1 == 6){
        offsetPosition.x = childSize;
        offsetPosition.y = childSize;
    } else if(checkOctant1 == 7){
        offsetPosition.z = childSize;
        offsetPosition.y = childSize;
        offsetPosition.x = childSize;
    }
    vec3 childBoxPosition = parentBox.position+offsetPosition;
    childBox.position = childBoxPosition;
    childBox.size = childSize;
    return childBox;
}


//Returns 0, 1, or 2. 0->has pushed, child found. 1->has pushed, no child found. 2->has not pushed, no child found.
uint push(rayStruct ray){
    uint checkOctant1 = octant1Stack[octant1StackQuantity-1];
    uint lastNodeIndex = nodeIndexStack[nodeIndexStackQuantity-1];
    uint lastNodeData = octree[lastNodeIndex];
    uint lastNodeValidMask8 = (lastNodeData >> 8) & 0x00FF;
    uint lastNodeLeafMask8 = lastNodeData & 0x00FF;
    if((lastNodeValidMask8 & (1 << 7-checkOctant1)) != 0){
        if((lastNodeLeafMask8 & (1 << 7-checkOctant1)) != 0){
            return 0; //continueMarching = false, successfullyPushed = true
        }
        uint lastNodeDataChildPointer = (lastNodeData >> 16) & 0xFFFF;
        uint newNodeAddress = calculateSiblingsBeforeThisOctant1(lastNodeData, checkOctant1) + lastNodeDataChildPointer + lastNodeIndex;
        if(nodeIndexStackQuantity < MAX_STACK_SIZE){
            nodeIndexStack[nodeIndexStackQuantity] = newNodeAddress;
            nodeIndexStackQuantity += 1;
        }
        
        uint newOctant1 = rayBoxOctant1(ray, boxStack[boxStackQuantity-1]);
        if(octant1StackQuantity < MAX_STACK_SIZE){
            octant1Stack[octant1StackQuantity] = newOctant1;
            octant1StackQuantity += 1;
        }
        boxStruct newBox = newBoxFromParentAndOctant1(boxStack[boxStackQuantity-1], octant1Stack[octant1StackQuantity-1]);
        if(boxStackQuantity < MAX_STACK_SIZE){
            boxStack[boxStackQuantity] = newBox;
            boxStackQuantity += 1;
        }
        return 1; //continueMarching = true, successfullyPushed = true
    }
    return 2; //continueMarching = true, successfullyPushed = false
}

uint advance(rayStruct ray) {
    boxStruct lastBox = boxStack[boxStackQuantity - 1u];
    boxStruct lastBoxParentBox = boxStack[boxStackQuantity - 2u];
    uint parentNodeIndex = nodeIndexStack[nodeIndexStackQuantity - 1u];
    uint parentNodeData = octree[parentNodeIndex];
    uint parentValidMask = (parentNodeData >> 8u) & 0x000000FF;
    uint parentLeafMask = parentNodeData & 0x000000FF;
    uint lastOctant1 = octant1Stack[octant1StackQuantity - 1u];
    bool leave = false;

    for (uint i = 0u; i < 4u; i++) {
        vec2 intersectData = checkBoxInverted(lastBox, ray);
        float tmax = intersectData.x;
        uint side = uint(intersectData.y);

        // Move to the new octant
        uint newOctant = moveOctant1(lastOctant1, side);

        if (newOctant == lastOctant1) {
            leave = true;
            break;
        }

        lastOctant1 = newOctant;
        lastBox = newBoxFromParentAndOctant1(lastBoxParentBox, lastOctant1);

        // Update the stacks
        boxStack[boxStackQuantity - 1u] = lastBox;
        octant1Stack[octant1StackQuantity - 1u] = lastOctant1;

        // Check if the ray has left the parent box
        vec2 parentIntersectData = checkBoxInverted(lastBoxParentBox, ray);
        float tmaxParent = parentIntersectData.x;

        if (tmax == tmaxParent) {
            leave = true;
        }
        
        if(leave == false){
            if ((parentValidMask & (1u << (7u - lastOctant1))) != 0u) {
                if ((parentLeafMask & (1u << (7u - lastOctant1))) != 0u) {
                    return 0u; // Found a leaf node
                } else {
                    return 1u; // Need to push
                }
            }
        }
        
    }

    if (leave) {
        return 2u; // Need to ascend to parent
    }

    return 1u; // Continue traversal
}




void splitPushCase(uint pushCase, inout bool continueMarching, inout bool successfullyPushed){
    if(pushCase == 0){
        continueMarching = false;
        successfullyPushed = true;
    } else if(pushCase == 1){
        continueMarching = true;
        successfullyPushed = true;
    } else if(pushCase == 2){
        continueMarching = true;
        successfullyPushed = false;
    }
}


void splitAdvanceCase(uint advanceCase, inout bool continueMarching, inout bool leaveChild){
    if(advanceCase == 0){
        leaveChild = false;
        continueMarching = false;
    } else if(advanceCase == 1){
        leaveChild = false;
        continueMarching = true;
    } else if(advanceCase == 2){
        leaveChild = true;
        continueMarching = true;
    }
}



vec2 marchRayThroughOctree(rayStruct ray, boxStruct boundingBox, vec3 uvInBoundingBox){
    bool continueMarching = true;
    bool pushed = false;
    bool leaveChild = false;

    uint firstChildOctant1 = calculateOctant1BasedOfEdgePos(uvInBoundingBox);
    boxStruct firstChildBox = newBoxFromParentAndOctant1(boundingBox, firstChildOctant1);
    if(boxStackQuantity < MAX_STACK_SIZE){
        boxStack[boxStackQuantity] = boundingBox;
        boxStackQuantity += 1;
    }
    if(boxStackQuantity < MAX_STACK_SIZE){
        boxStack[boxStackQuantity] = firstChildBox;
        boxStackQuantity += 1;
    }
    if(nodeIndexStackQuantity < MAX_STACK_SIZE){
        nodeIndexStack[nodeIndexStackQuantity] = 0;
        nodeIndexStackQuantity += 1;
    }
    if(octant1StackQuantity < MAX_STACK_SIZE){
        octant1Stack[octant1StackQuantity] = firstChildOctant1;
        octant1StackQuantity += 1;
    }
    
    int raySteps = 0;
    for(int rayStep = 0; rayStep < MAX_RAY_STEPS; rayStep++){
        if(!continueMarching){
            break;
        }
        uint pushCase = push(ray);
        if(pushCase == 0){
            continueMarching = false;
            pushed = true;
        } else if(pushCase == 1){
            continueMarching = true;
            pushed = true;
        } else if(pushCase == 2){
            continueMarching = true;
            pushed = false;
        }
        
        if(!pushed){
            uint advanceCase = advance(ray);
            splitAdvanceCase(advanceCase, continueMarching, leaveChild);
            while(leaveChild == true){
                if(boxStackQuantity > 0){
                    boxStackQuantity -= 1;
                    boxStruct emptyBox;
                    emptyBox.position = vec3(0.0, 0.0, 0.0);
                    emptyBox.size = 0.0;
                    boxStack[boxStackQuantity] = emptyBox;
                }
                if(nodeIndexStackQuantity > 0){
                    nodeIndexStackQuantity -= 1;
                    nodeIndexStack[nodeIndexStackQuantity] = 0;
                }
                if(octant1StackQuantity > 0){
                    octant1StackQuantity -= 1;
                    octant1Stack[octant1StackQuantity] = 0;
                }
                advanceCase = advance(ray);
                if(advanceCase == 0){
                    leaveChild = false;
                    continueMarching = false;
                } else if(advanceCase == 1){
                    leaveChild = false;
                    continueMarching = true;
                } else if(advanceCase == 2){
                    leaveChild = true;
                    continueMarching = true;
                }
                if(octant1StackQuantity < 1){
                    return vec2(0, raySteps); //Blue
                }
            }
            
        } 
        raySteps += 1;
    }

    
    boxStruct foundBox = boxStack[boxStackQuantity-1];
    vec2 foundBoxIntersect = checkBox(foundBox, ray);
    float t = foundBoxIntersect[0];
    if(t < 0){
        return vec2(0, raySteps); //Green
    }        
    return vec2(t, raySteps);    
    
}


void createBox(vec3 position, float size, int index){
    boxStruct newBox;
    newBox.position = position;
    newBox.size = size;
    boxes[index] = newBox;
}


vec4 checkAllBoxes(rayStruct ray){
    vec2 intersectData;
    vec2 temporaryIntersectData;
    int intersectCount = 0;
    intersectData[0] = 1e30;
    intersectData[1] = -1e30;
    for(int i = 0; i < boxes.length(); i++){
        temporaryIntersectData = checkBoxInverted(boxes[i], ray);
        if(boxes[i].size != 0.0){
            if(temporaryIntersectData[0] > 0 && temporaryIntersectData[0] < intersectData[0]){
                intersectData = temporaryIntersectData;
                intersectCount += 1;
            }
        }
        if(temporaryIntersectData[0] > 0){
            intersectCount += 1;
        }
        
    }
    //if(intersectData[0] >= 1e29){ intersectData[0] = -1; }
    return vec4(intersectData[0], intersectData[1], intersectCount, 0.0);
}


vec4 castRay(){

    vec4 color;
    rayStruct ray;

    ray.direction = rayStartDirection(TexCoords, resolution, cameraPos, cameraDir, 90);
    ray.origin = cameraPos;

    boxStruct octreeBoundingBox;
    octreeBoundingBox.position = vec3(20.0, 0.0, 0.0);
    octreeBoundingBox.size = 20.0;
    vec2 octreeBoundingBoxIntersectData = checkBox(octreeBoundingBox, ray);
    float tBB = octreeBoundingBoxIntersectData[0];
    vec3 intersectionPoint = ray.origin+ray.direction*tBB;

    vec3 uvInBoundingBox;
    if (tBB < 0.0) { // Camera is inside the bounding box
        uvInBoundingBox = abs((octreeBoundingBox.position - ray.origin)) / octreeBoundingBox.size;
    } else {
        vec3 intersectionPoint = ray.origin + ray.direction * tBB;
        uvInBoundingBox = abs((octreeBoundingBox.position - intersectionPoint)) / octreeBoundingBox.size;
    }

    
    vec2 octreeMarchData;

    if(max(max(uvInBoundingBox.x, uvInBoundingBox.y), uvInBoundingBox.z) <= 1 && min(min(uvInBoundingBox.x, uvInBoundingBox.y), uvInBoundingBox.z) >= 0){
        
        octreeMarchData = marchRayThroughOctree(ray, octreeBoundingBox, uvInBoundingBox);
    }
    
    if(octreeMarchData[1]/MAX_RAY_STEPS == 1){
        return vec4(0.9, 0.0, 0.0, 1.0);
    }
    //return vec4(octreeMarchData.x/20, octreeMarchData.x/20, octreeMarchData.x/20, 1.0);
    vec2 temp = octreeBoundingBoxIntersectData;
    
    //return vec4(0.0, 0.0, 0.0, 0.0);
    return vec4(octreeMarchData[1]/MAX_RAY_STEPS, octreeMarchData[1]/MAX_RAY_STEPS, octreeMarchData[1]/MAX_RAY_STEPS, 1.0);
    if(octreeMarchData[1] > 1){
        //return vec4(octreeMarchData[0]/10, octreeMarchData[0]/10, octreeMarchData[0]/10, 1.0);
    }
}


void main(){
    FragColor = castRay();
}
