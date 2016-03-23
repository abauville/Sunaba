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
#include <float.h>







//============================================================================//
//============================================================================//
//                                                                            //
//                             COMPILE SHADER FUNCTION                        //
//                                                                            //
//============================================================================//
//============================================================================//

void compileShaders(GLuint *ShaderProgram, const char* pVSFileName, const char* pFSFileName, const char* pGSFileName, bool useGS)
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
	char *vs, *fs, *gs;
	vs = readFile(pVSFileName);
	fs = readFile(pFSFileName);
	if (useGS)
		gs = readFile(pGSFileName);




	// =========================================================================
	//                                 SHADERS
	// =========================================================================

	// Init shaders
	// =======================================
	GLuint VertexShader = glCreateShader(GL_VERTEX_SHADER);
	GLuint FragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

	GLuint GeometryShader;
	if (useGS) {
		GeometryShader= glCreateShader(GL_GEOMETRY_SHADER);
	}




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
	if (useGS) {
		p[0] = gs;    Lengths[0] = strlen(gs);
		glShaderSource(GeometryShader,1,p,Lengths);
		glCompileShader(GeometryShader);
	}


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
	// GeometryShader
	if (useGS) {
		glGetShaderiv(GeometryShader, GL_COMPILE_STATUS, &success);
		if (!success) {
			GLchar InfoLog[1024];
			glGetShaderInfoLog(GeometryShader, sizeof(InfoLog), NULL, InfoLog);
			fprintf(stderr, "Error compiling shader type %d: '%s'\n", GL_FRAGMENT_SHADER, InfoLog);
		}
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
	if (useGS) {
		glAttachShader(*ShaderProgram, GeometryShader);
	}
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
//                             MIN, MAX and friends                           //
//                                                                            //
//============================================================================//
//============================================================================//
double min(double* List, int length)
{
	int i;
	double MIN = DBL_MAX; // largest number
	for (i = 0; i < length; ++i) {
		if (List[i]<MIN)
			MIN = (List[i]);
	}
	return MIN;
}
double max(double* List, int length)
{
	int i;
	double MAX = DBL_MIN;
	for (i = 0; i < length; ++i) {
		if (List[i]>MAX) {
			MAX = (List[i]);
		}
	}

	return MAX;
}

float minf(float* List, int length)
{
	int i;
	float MIN = FLT_MAX; // largest number
	for (i = 0; i < length; ++i) {
		if (List[i]<MIN)
			MIN = (List[i]);
	}
	return MIN;
}
float maxf(float* List, int length)
{
	int i;
	float MAX = FLT_MIN;
	for (i = 0; i < length; ++i) {
		if (List[i]>MAX)
			MAX = (List[i]);
	}
	return MAX;
}

double absmin(double* List, int length)
{
	int i;
	double MIN = DBL_MAX; // largest number
	for (i = 0; i < length; ++i) {
		if (fabs(List[i])<MIN)
			MIN = fabs(List[i]);
	}
	return MIN;
}
double absmax(double* List, int length)
{
	int i;
	double MAX = 0;
	for (i = 0; i < length; ++i) {
		if (fabs(List[i])>MAX)
			MAX = fabs(List[i]);
	}
	return MAX;
}






//============================================================================//
//============================================================================//
//                                                                            //
//                                BITMAP READER                               //
//                                                                            //
// This source code is written by Mikito Furuichi since 2013/03/24  		  //
// mfuruic@gmail.com                                                          //
//============================================================================//
//============================================================================//

void ReadBmp(char *filename, img *imgp) {
  int i,j;
  int Real_width;
  FILE *Bmp_Fp=fopen(filename,"rb");
  unsigned char *Bmp_Data;

  if(Bmp_Fp==NULL){
    fprintf(stderr,"Error: file %s couldn\'t open for read!.\n",filename);
    exit(1);
  }

  fread(Bmp_headbuf,sizeof(unsigned char),HEADERSIZE,Bmp_Fp);

  memcpy(&Bmp_type,Bmp_headbuf,sizeof(Bmp_type));
  if (strncmp(Bmp_type,"BM",2)!=0) {
    fprintf(stderr,"Error: %s is not a bmp file.\n",filename);
    exit(1);
  }

  memcpy(&imgp->width,Bmp_headbuf+18,sizeof(Bmp_width));
  memcpy(&imgp->height,Bmp_headbuf+22,sizeof(Bmp_height));
  memcpy(&Bmp_color,Bmp_headbuf+28,sizeof(Bmp_color));
  if (Bmp_color!=24) {
    fprintf(stderr,"Error: Bmp_color = %d is not implemented in this program.\n",Bmp_color);
    exit(1);
  }

  if (imgp->width > MAXWIDTH) {
    fprintf(stderr,"Error: Bmp_width = %ld > %d = MAXWIDTH!\n",Bmp_width,MAXWIDTH);
    exit(1);
  }

  if (imgp->height > MAXHEIGHT) {
    fprintf(stderr,"Error: Bmp_height = %ld > %d = MAXHEIGHT!\n",Bmp_height,MAXHEIGHT);
    exit(1);
  }

  Real_width = imgp->width*3 + imgp->width%4;

 if((Bmp_Data = (unsigned char *)calloc(Real_width,sizeof(unsigned char)))==NULL) {
   fprintf(stderr,"Error: Memory allocation failed for Bmp_Data!\n");
   exit(1);
 }

  for(i=0;i<imgp->height;i++) {
    fread(Bmp_Data,1,Real_width,Bmp_Fp);
    for (j=0;j<imgp->width;j++) {
      imgp->data[imgp->height-i-1][j].b = Bmp_Data[j*3];
      imgp->data[imgp->height-i-1][j].g = Bmp_Data[j*3+1];
      imgp->data[imgp->height-i-1][j].r = Bmp_Data[j*3+2];
    }
  }

  free(Bmp_Data);

  fclose(Bmp_Fp);
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



