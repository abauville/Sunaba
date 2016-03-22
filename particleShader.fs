#version 330

in vec3 fColor;
out vec4 out_Color;

void main() {
    
    out_Color = vec4(fColor,1.0);
}