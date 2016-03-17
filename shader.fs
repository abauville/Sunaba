#version 330 core

in vec3 color;
out vec4 out_Color;

void main() {
    
    
    //set every drawn pixel to white
    out_Color = vec4(color, 1.0);
}