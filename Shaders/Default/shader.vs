#version 330 core


in vec2 in_Vertex;
in vec2 in_TexCoord;
//in float U;

//out vec3 color;
out vec2 TexCoord;



uniform mat4 transform;


void main() {
    

    TexCoord = in_TexCoord;
    gl_Position = transform * vec4(in_Vertex, 0.0, 1);

    //color = texture(tex, in_TexCoord) * vec4(1.0,1.0,1.0,1.0);
   // color = vec4(1.0,1.0,1.0,1.0);

   // color = vec3(1.0,1.0,1.0);
    

}