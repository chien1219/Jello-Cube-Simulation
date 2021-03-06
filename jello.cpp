/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

  Your name: Li-Che (Richard) Chien

*/

#include "jello.h"
#include "showCube.h"
#include "input.h"
#include "physics.h"
#include "MyFunc.h"
#include "Texture\Cbmp.h"

bool firsttime=true;

double sw, sh;
double tw, th;

GLuint g_uiTextureId[g_ciTexNum] = {0};

// camera parameters
double Theta = pi / 6;
double Phi = pi / 6;
double R = 6;

// mouse control
int g_iMenuId;
int g_vMousePos[2];
int g_iLeftMouseButton,g_iMiddleMouseButton,g_iRightMouseButton;
int DownPos[2];
int UpPos[2];
bool isMouseForce;

// number of images saved to disk so far
int sprite=0;

double JelloFixedR;
int VerticesState[8][8][8];

int shear = 0, bend = 0, structural = 1, pause = 0, viewingMode = 1, saveScreenToFile = 0, showTexture = 0, EXforce = 1;

struct world jello;

int windowWidth, windowHeight;

void TextureInit()
{
	CBmp TextureImage[g_ciTexNum];
	TextureImage[0].load("Texture/dice0.bmp");
	TextureImage[1].load("Texture/dice1.bmp");
	TextureImage[2].load("Texture/dice2.bmp");
	TextureImage[3].load("Texture/dice3.bmp");
	TextureImage[4].load("Texture/dice4.bmp");
	TextureImage[5].load("Texture/dice5.bmp");

	for (int iI = 0; iI < g_ciTexNum; iI++)
	{
		glEnable(GL_TEXTURE_2D);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glGenTextures(1, &g_uiTextureId[iI]);
		glBindTexture(GL_TEXTURE_2D, g_uiTextureId[iI]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, TextureImage[iI].w, TextureImage[iI].h, 0, GL_RGB, GL_UNSIGNED_BYTE, TextureImage[iI].rgb);
		//gluBuild2DMipmaps(GL_TEXTURE_2D, 3 , TextureImage[iI].w, TextureImage[iI].h,GL_RGB ,GL_UNSIGNED_BYTE, TextureImage[iI].rgb);
		glBindTexture(GL_TEXTURE_2D, 0);
		glDisable(GL_TEXTURE_2D);
	}
}

void myinit()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(90.0,1.0,0.01,1000.0);
	
	// set background color to grey
	glClearColor(0.5, 0.5, 0.5, 0.0);
	
	glCullFace(GL_BACK);
	glEnable(GL_CULL_FACE);
	
	glShadeModel(GL_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_LINE_SMOOTH);

	 // Texture
    if(firsttime)
    {
		Vector L;
		pDIFFERENCE(jello.p[0][0][0], jello.p[0][0][1], L);
		JelloFixedR = ComputeLength(L);

		firsttime = false;

		TextureInit();
    }
  return;
}

void reshape(int w, int h)
{
	// Prevent a divide by zero, when h is zero.
	// You can't make a window of zero height.
	if(h == 0)
	  h = 1;
	
	glViewport(0, 0, w, h);
	
	// Reset the coordinate system before modifying
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	// Set the perspective
	double aspectRatio = 1.0 * w / h;
	gluPerspective(60.0f, aspectRatio, 0.01f, 1000.0f);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	windowWidth = w;
	windowHeight = h;
	
	glutPostRedisplay();
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	// camera parameters are Phi, Theta, R
	gluLookAt(R * cos(Phi) * cos (Theta), R * sin(Phi) * cos (Theta), R * sin (Theta),
	        0.0,0.0,0.0, 0.0,0.0,1.0);
	
	
	/* Lighting */
	/* You are encouraged to change lighting parameters or make improvements/modifications
	   to the lighting model .
	   This way, you will personalize your assignment and your assignment will stick out.
	*/
	
	// global ambient light
	GLfloat aGa[] = { 0.0, 0.0, 0.0, 0.0 };
	
	// light 's ambient, diffuse, specular
	GLfloat lKa0[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd0[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat lKs0[] = { 1.0, 1.0, 1.0, 1.0 };
	
	GLfloat lKa1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd1[] = { 1.0, 0.0, 0.0, 1.0 };
	GLfloat lKs1[] = { 1.0, 0.0, 0.0, 1.0 };
	
	GLfloat lKa2[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd2[] = { 1.0, 1.0, 0.0, 1.0 };
	GLfloat lKs2[] = { 1.0, 1.0, 0.0, 1.0 };
	
	GLfloat lKa3[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd3[] = { 0.0, 1.0, 1.0, 1.0 };
	GLfloat lKs3[] = { 0.0, 1.0, 1.0, 1.0 };
	
	GLfloat lKa4[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd4[] = { 0.0, 0.0, 1.0, 1.0 };
	GLfloat lKs4[] = { 0.0, 0.0, 1.0, 1.0 };
	
	GLfloat lKa5[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd5[] = { 1.0, 0.0, 1.0, 1.0 };
	GLfloat lKs5[] = { 1.0, 0.0, 1.0, 1.0 };
	
	GLfloat lKa6[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd6[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat lKs6[] = { 1.0, 1.0, 1.0, 1.0 };
	
	GLfloat lKa7[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd7[] = { 0.0, 1.0, 1.0, 1.0 };
	GLfloat lKs7[] = { 0.0, 1.0, 1.0, 1.0 };
	
	// light positions and directions
	GLfloat lP0[] = { -1.999, -1.999, -1.999, 1.0 };
	GLfloat lP1[] = { 1.999, -1.999, -1.999, 1.0 };
	GLfloat lP2[] = { 1.999, 1.999, -1.999, 1.0 };
	GLfloat lP3[] = { -1.999, 1.999, -1.999, 1.0 };
	GLfloat lP4[] = { -1.999, -1.999, 1.999, 1.0 };
	GLfloat lP5[] = { 1.999, -1.999, 1.999, 1.0 };
	GLfloat lP6[] = { 1.999, 1.999, 1.999, 1.0 };
	GLfloat lP7[] = { -1.999, 1.999, 1.999, 1.0 };
	
	// jelly material color
	GLfloat mKa[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat mKd[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat mKs[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mKe[] = { 0.3, 0.3, 0.3, 1.0 };
	
	/* set up lighting */
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, aGa);
	glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
	
	// set up cube color
	glMaterialfv(GL_FRONT, GL_AMBIENT, mKa);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mKd);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mKs);
	glMaterialfv(GL_FRONT, GL_EMISSION, mKe);
	glMaterialf(GL_FRONT, GL_SHININESS, 120);
	
	// macro to set up light i
	#define LIGHTSETUP(i)\
	glLightfv(GL_LIGHT##i, GL_POSITION, lP##i);\
	glLightfv(GL_LIGHT##i, GL_AMBIENT, lKa##i);\
	glLightfv(GL_LIGHT##i, GL_DIFFUSE, lKd##i);\
	glLightfv(GL_LIGHT##i, GL_SPECULAR, lKs##i);\
	glEnable(GL_LIGHT##i)
	
	LIGHTSETUP (0);
	LIGHTSETUP (1);
	LIGHTSETUP (2);
	LIGHTSETUP (3);
	LIGHTSETUP (4);
	LIGHTSETUP (5);
	LIGHTSETUP (6);
	LIGHTSETUP (7);
	
	// enable lighting
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
	
	// show the cube
	showCube(&jello);
	
	//glDisable(GL_LIGHTING);
	
	// show the bounding box
	showBoundingBox();
	
	glutSwapBuffers();
}

void doIdle()
{
  char s[20]="picxxxx.ppm";
  int i;

  // save screen to file
  s[3] = 48 + (sprite / 1000);
  s[4] = 48 + (sprite % 1000) / 100;
  s[5] = 48 + (sprite % 100 ) / 10;
  s[6] = 48 + sprite % 10;

  if (saveScreenToFile==1)
  {
    saveScreenshot(windowWidth, windowHeight, s);
    //saveScreenToFile=0; // save only once, change this if you want continuos image generation (i.e. animation)
    sprite++;
  }

  if (sprite >= 300) // allow only 300 snapshots
  {
    exit(0);
  }

  if (pause == 0)
  {
      // insert code which appropriately performs one step of the cube simulation:
      if (strcmp(jello.integrator, "Euler") == 0)
		  Euler(&jello);
      else
		  RK4(&jello);
  }

  glutPostRedisplay();
}

void InitAllPointsSide(struct world * jello)
{
	int i,j,k;
	
	if(jello->incPlanePresent)
	{
		for (i=0; i<=7; i++)
			for (j=0; j<=7; j++)
				for (k=0; k<=7; k++)
					VerticesState[i][j][k]=jello->a*jello->p[i][j][k].x+ jello->b*jello->p[i][j][k].y+ jello->c*jello->p[i][j][k].z+ jello->d;
	}
}

int main (int argc, char ** argv)
{
  if (argc<2)
  {
    printf ("Oops! You didn't say the jello world file!\n");
    printf ("Usage: %s [worldfile]\n", argv[0]);
    exit(0);
  }

  readWorld(argv[1],&jello);
  InitAllPointsSide(&jello);


  glutInit(&argc,argv);

  /* double buffered window, use depth testing, 640x480 */
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

  windowWidth = 640;
  windowHeight = 480;
  glutInitWindowSize (windowWidth, windowHeight);
  glutInitWindowPosition (0,0);
  glutCreateWindow ("Jello cube");

  /* tells glut to use a particular display function to redraw */
  glutDisplayFunc(display);

  /* replace with any animate code */
  glutIdleFunc(doIdle);

  /* callback for mouse drags */
  glutMotionFunc(mouseMotionDrag);

  /* callback for window size changes */
  glutReshapeFunc(reshape);

  /* callback for mouse movement */
  glutPassiveMotionFunc(mouseMotion);

  /* callback for mouse button changes */
  glutMouseFunc(mouseButton);

  /* register for keyboard events */
  glutKeyboardFunc(keyboardFunc);

  /* do initialization */
  myinit();

  /* forever sink in the black hole */
  glutMainLoop();

  return(0);
}

