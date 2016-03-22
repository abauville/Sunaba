#version 330

in vec2 PartVertex;
in float PartData;
uniform mat4 transform;
out vec3 vColor;


void main() {
    

    gl_Position = transform * vec4(PartVertex, 0.0, 1);
    //gl_Position = vec4(PartVertex, 0.0, 1);
    
    if (PartData == 0) {
        vColor = vec3(0.3,1.0,0.5);
    }
    else if (PartData == 1) {
        vColor = vec3(1.0,0.6,0.3);
    }
    else if (PartData == 2) {
        vColor = vec3(0.0,1.0,0.0);
    }
    else if (PartData == 3) {
        vColor = vec3(0.0,0.0,1.0);
    }
    else if (PartData == 4) {
        vColor = vec3(1.0,1.0,0.0);
    }
    else if (PartData == 5) {
        vColor = vec3(0.0,1.0,1.0);
    }
    else {
        vColor = vec3(1.0,1.0,1.0);
    }
    
    

}