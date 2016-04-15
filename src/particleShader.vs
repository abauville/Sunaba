#version 330

in vec2 PartVertex;
in vec3 PartMeshVertex;
in float PartData;
in float PartPassiveData;
uniform mat4 transform;
out vec3 vColor;
uniform float size;


void main() {
    //gl_Position = transform * vec4(PartVertex, 0.0, 1);
    vec4 Pos = ( vec4(PartVertex, 0.0, 1) + vec4(PartMeshVertex, 0) );
    gl_Position = transform*Pos;
    
    if (PartData == 0) {
        vColor = vec3(0.0,0.0,0.0);
    }
    else if (PartData == 1) {
        vColor = vec3(0.3,1.0,0.5);
    }
    else if (PartData == 2) {
        vColor = vec3(1.0,0.5,0.3);
    }
    else if (PartData == 3) {
        vColor = vec3(0.3,0.5,1.0);
    }
    else if (PartData == 4) {
        vColor = vec3(0.7,1.0,0.1);
    }
    else if (PartData == 5) {
        vColor = vec3(0.0,1.0,1.0);
    }
    else {
        vColor = vec3(1.0,1.0,1.0);
    }
    
    
    if (PartPassiveData<0.5) {
        vColor = 0.6*vColor;
    }
    

}