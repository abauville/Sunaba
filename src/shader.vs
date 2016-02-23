#version 330 core


in vec2 in_Vertex;
in float U;

out vec3 color;
uniform mat4 transform;

void main() {
    // Final Position of the Vertices
    gl_Position = transform * vec4(in_Vertex, 0.0, 1);

    // Compute values of color according to the solution
    float ca1, ca2, cb1, cb2, cc1, cc2, cd1, cd2; // limit of caxis for each color
    float R, G, B;
    //ca1 = 0.5;      ca2 = 1;
    //cb1 = 1;        cb2 = 2.5;
    //cc1 = -0.5;     cc2 = -1;
    
    
    ca1 = 0.0;      ca2 =  .5;
    cb1 = ca2;    cb2 =  2*ca2;
    cc1 = -ca1;   cc2 = -ca2;
    
    
    R = (U-ca1)/(ca2-ca1);
    G = (U-cb1)/(cb2-cb1);
    B = (U-cc1)/(cc2-cc1);
    
    if (G<0.0) {
        G = (U+cb1)/(-(cb2-cb1));
    }
    
    R = U;
    
    if (R>1.0){
        R = 1.0;
    }
    if (G>1.0){
        G = 1.0;
    }
    if (B>1.0){
        B = 1.0;
    }
    
    
    if (R<0.0){
        R = 0.0;
    }
    if (G<0.0){
        G = 0.0;
    }
    if (B<0.0){
        B = 0.0;
    }
    
    
    // Send color to the fragment shader
    color.xyz = vec3(R, G, B);
        
    

}