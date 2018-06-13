#version 330 core

layout(points) in;
layout(triangle_strip, max_vertices = 5) out;

in vec4 vColor[]; // Output from vertex shader for each vertex

out vec3 fColor; // Output to fragment shader

uniform float size;

void main() {
    //float size = 0.1;
    fColor = vColor[0];
    
    gl_Position = gl_in[0].gl_Position + vec4(-size, -size, 0.1, 0.0);
    EmitVertex();
    
    gl_Position = gl_in[0].gl_Position + vec4(0.0, size, 0.1, 0.0);
    EmitVertex();
    
    gl_Position = gl_in[0].gl_Position + vec4(0.0,0.0, -0.1, 0.0); // center node
    EmitVertex();
    
    gl_Position = gl_in[0].gl_Position + vec4(size, -size, 0.1, 0.0);
    EmitVertex();
    
    gl_Position = gl_in[0].gl_Position + vec4(-size, -size, 0.1, 0.0);
    EmitVertex();

    
    EndPrimitive();
}