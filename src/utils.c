//
//  utils.c
//
//
//  Created by Arthur Bauville on 1/26/16.
//
//

#include "stokes.h"
#include <string.h>
#include <stdlib.h>







//============================================================================//
//============================================================================//
//                                                                            //
//                             COMPILE SHADER FUNCTION                        //
//                                                                            //
//============================================================================//
//============================================================================//

void compileShaders(GLuint *ShaderProgram, const char* pVSFileName, const char* pFSFileName)
{
    // =========================================================================
    //                                  INIT
    // =========================================================================

    // Initialize variables
    // =======================================
    const GLchar* p[1];
    GLint Lengths[1];
    GLchar ErrorLog[1024] = { 0 };
    GLint success;


    // Read shader file and put it in a string
    // =======================================
    char *vs, *fs;
    vs = readFile(pVSFileName);
    fs = readFile(pFSFileName);





    // =========================================================================
    //                                 SHADERS
    // =========================================================================

    // Init shaders
    // =======================================
    GLuint VertexShader = glCreateShader(GL_VERTEX_SHADER);
    GLuint FragmentShader = glCreateShader(GL_FRAGMENT_SHADER);


    // Compile shaders
    // =======================================
    // VertexShader
    p[0] = vs;    Lengths[0] = strlen(vs);
    glShaderSource(VertexShader,1,p,Lengths);
    glCompileShader(VertexShader);
    // FragmentShader
    p[0] = fs;    Lengths[0] = strlen(fs);
    glShaderSource(FragmentShader,1,p,Lengths);
    glCompileShader(FragmentShader);


    // Check for compilation errors
    // =======================================
    // VertexShader
    glGetShaderiv(VertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        GLchar InfoLog[1024];
        glGetShaderInfoLog(VertexShader, sizeof(InfoLog), NULL, InfoLog);
        fprintf(stderr, "Error compiling shader type %d: '%s'\n", GL_VERTEX_SHADER, InfoLog);
    }
    // FragmentShader
    glGetShaderiv(FragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        GLchar InfoLog[1024];
        glGetShaderInfoLog(FragmentShader, sizeof(InfoLog), NULL, InfoLog);
        fprintf(stderr, "Error compiling shader type %d: '%s'\n", GL_FRAGMENT_SHADER, InfoLog);
    }





    // =========================================================================
    //                                 PROGRAM
    // =========================================================================

    // Create program and attach shader
    // =======================================
    *ShaderProgram = glCreateProgram();


    // Attach shaders
    // =======================================
    glAttachShader(*ShaderProgram, VertexShader);
    glAttachShader(*ShaderProgram, FragmentShader);


    // Link the program
    // =======================================
    /// Link the shader program to the rest of the program
    glLinkProgram(*ShaderProgram);


    // Check for linking errors
    // =======================================
    glGetProgramiv(*ShaderProgram, GL_LINK_STATUS, &success);
    if (success == 0) {
        glGetProgramInfoLog(*ShaderProgram, sizeof(ErrorLog), NULL, ErrorLog);
        fprintf(stderr, "Error linking shader program: '%s'\n", ErrorLog);
        exit(1);
    }

    // Validate the program
    // =======================================
    // (can it execute given the current pipeline?)
    glValidateProgram(*ShaderProgram);
    glGetProgramiv(*ShaderProgram, GL_VALIDATE_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(*ShaderProgram, sizeof(ErrorLog), NULL, ErrorLog);
        fprintf(stderr, "Invalid shader program: %s\n", ErrorLog);
        exit(1);
    }


    // =========================================================================
    //                                 FREE
    // =========================================================================

    // Detach shaders and free code source strings
    // =======================================
    glDetachShader(*ShaderProgram, VertexShader);
    glDetachShader(*ShaderProgram, FragmentShader);
    glDetachShader(*ShaderProgram,VertexShader);
    glDetachShader(*ShaderProgram,FragmentShader);
    free(vs);
    free(fs);



}



//============================================================================//
//============================================================================//
//                                                                            //
//                               READ FILE FUNCTION                           //
//                                                                            //
//============================================================================//
//============================================================================//

char* readFile(const char* fileName)
{
    GLchar* myFileString;
    FILE* fp = NULL; // file pointer
    int size;


    fp = fopen(fileName,"r");
    if (fp == NULL)
    {
        fprintf(stderr, "error:  %s: file not found.\n", fileName);
        exit(1);
    }

    // Get the file's size
    fseek(fp, 0L, SEEK_END);
    size = ftell(fp);
    fseek(fp, 0L, SEEK_SET);

    // Allocate memory for the file string
    myFileString = (GLchar*) malloc(size * sizeof(GLchar));

    fread (myFileString, 1, size, fp);
    myFileString[size] = 0; // NULL terminator

    fclose(fp);
    return myFileString;

}




//============================================================================//
//============================================================================//
//                                                                            //
//                                 LIST PRINTING                              //
//                                                                            //
//============================================================================//
//============================================================================//

void printListi(int* list, int length) {
	int i;
	for(i=0;i<length;i++){
		printf("%i  ", list[i]);
	}
	printf("\n");
}


void printListd(double* list,int length) {
	int i;
	for(i=0;i<length;i++){
		printf("%.1f  ", list[i]);
	}
	printf("\n");
}





//============================================================================//
//============================================================================//
//                                                                            //
//                                    FIND                                    //
//                                                                            //
//============================================================================//
//============================================================================//
int findi(int* list, int length, int value)
{
	int i = 0;
	while ( (list[i] != value) && (i<length) )
	{
		i++;
	}
	if (i==length) {
		i = -1;
	}
	return i;
}





//============================================================================//
//============================================================================//
//                                                                            //
//                                   GRAVEYARD                                //
//                                                                            //
//============================================================================//
//============================================================================//


// Validate the program
// =======================================
/*
 // Validate the program (can it execute given the current pipeline?)
 glValidateProgram(*ShaderProgram);
 glGetProgramiv(*ShaderProgram, GL_VALIDATE_STATUS, &success);
 if (!success) {
 glGetProgramInfoLog(*ShaderProgram, sizeof(ErrorLog), NULL, ErrorLog);
 int ret;
 ret = strcmp(ErrorLog,"Validation Failed: No vertex array object bound.\n");
 if (ret==0) {
 fprintf(stderr,"Warning: Validation Failed: No vertex array object bound. (Strange Mac OS X bug)\n");
 }
 else {
 fprintf(stderr, "Invalid shader program: %s\n", ErrorLog);
 exit(1);
 }
 }
 */

/// Load the shader in the pipeline
//glUseProgram(ShaderProgram);



