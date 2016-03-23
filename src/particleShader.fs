#version 330

in vec3 vColor;
out vec4 out_Color;

void main() {
    
    out_Color = vec4(vColor,1.0);
}