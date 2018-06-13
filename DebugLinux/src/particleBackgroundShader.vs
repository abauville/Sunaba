#version 330 core


in vec2 in_Vertex;

uniform mat4 transform;
out vec3 color;

void main() {
    
    gl_Position = transform * vec4(in_Vertex, +0.0001, 1.0);


    color = vec3(0.0,0.0,0.0);
    

}