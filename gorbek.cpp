#include <cmath>
#include <vector>

#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include "bevgrafmath2017.h"

//	Size of window
GLsizei winHeight = 800, winWidth = 1000;

std::vector<vec2> points = { {66, 786}, {44, 653}, {285, 776}, {578, 738},{703, 719},
							{834, 445},{466, 475},{190, 497.5},{42, 349},
							{171, 136},{507, 182} };

vec2 erintoPont = {0,0};

std::vector<vec2> originalPoints = points;

//	Hermite-ív paraméter tartománya
GLfloat start = -2;
GLfloat middle = 0.5;
GLfloat end = 1.5;


//////////////////////////////  Hermite-ív	///////////////////////////////////
mat24 G = {points[1], points[2], points[3], points[0]-points[1]};

GLfloat t1 = start;
GLfloat t2 = middle;
GLfloat t3 = end;

vec4 T1 = {t1*t1*t1, t1*t1, t1, 1};
vec4 T2 = {t2*t2*t2, t2*t2, t2, 1};
vec4 T3 = {t3*t3*t3, t3*t3, t3, 1};
vec4 T4 = {3*t1*t1, 2*t1, 1, 0};

mat4 M = inverse({T1,T2,T3,T4, true});

vec4 TT = {3*t3*t3, 2*t3, 1, 0};

///////////////////////////////////////////////////////////////////////////////

//	de Casteljau parameter
GLfloat u = 0.5;

// set up pick radius for detecting movement of a control point
GLint pickRadius = 4;

//	aktuális pont kiválasztása
GLint dragged = -1;

//	show/hide control polygon
bool displayPoligon = true;

void calculatePoints() {
	points[4] = points[3] + (G * M * TT)/3;
	erintoPont = points[3] + (G * M * TT);
	points[7] = points[6] +  0.75 * (points[6] - points[5]);
}

void displayControlPolygon() {
	calculatePoints();

	glLineWidth(1.0);
	glColor3f(0.5, 0.5, 0.5);

	glBegin(GL_LINE_STRIP);
	for (int i = 0; i < 11; i++)
	{
	  glVertex2f(points[i].x,points[i].y);
	}
	glEnd();

	glBegin(GL_LINE_STRIP);
	glVertex2f(erintoPont.x, erintoPont.y);
	glVertex2f(points[4].x,points[4].y);
	glEnd();
}

void displayControlPoints() {
	calculatePoints();

	glPointSize(8);
	glBegin(GL_POINTS);
	for (int i = 0; i < 11; i++)
	{
		glColor3f(1.0, 0.0, 1.0);
		if(i == 4 || i == 7)
			glColor3f(0.0, 1.0, 1.0);
		glVertex2f(points[i].x,points[i].y);
	}
	glColor3f(0.0, 1.0, 1.0);
	glVertex2f(erintoPont.x, erintoPont.y);
	glEnd();

}

void hermite() {
	calculatePoints();
	G = {points[1], points[2], points[3], points[0]-points[1]};

	glLineWidth(3.0);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINE_STRIP);
	for (float t = start; t <= end; t+=0.001) {
		vec4 T = {t*t*t, t*t, t, 1};
		vec2 curvePoint = G*M*T;
		glVertex2f(curvePoint.x, curvePoint.y);
	}
	glEnd();
}

void bernstein() {
	calculatePoints();

	glLineWidth(3.0);
	glColor3f(0.0 , 1.0, 0.0);
	glBegin(GL_LINE_STRIP);

	for (float t = 0; t <= 1; t+=0.001) {
		vec2 curvePoint =
		(-(pow(t,3)) + 3*pow(t,2) -3*t + 1) * points[3] +
		(3*pow(t,3) - 6*pow(t,2) + 3*t) * points[4] +
		(-3*pow(t,3) + 3*pow(t,2)) * points[5] +
		pow(t,3) * points[6];

		glVertex2f(curvePoint.x, curvePoint.y);
	}
	glEnd();
}

vec2 oszto_pont(vec2 a, vec2 b, GLfloat t) {
	GLfloat s = (1-t);
	return {a.x*s + b.x*t, a.y*s + b.y*t };
}

vec2 calculateCurvePoint(std::vector<vec2> pontok, GLfloat uu) {

	vec2 b0[4];
	b0[0] = oszto_pont(pontok[0], pontok[1], uu);
	b0[1] = oszto_pont(pontok[1], pontok[2], uu);
	b0[2] = oszto_pont(pontok[2], pontok[3], uu);
	b0[3] = oszto_pont(pontok[3], pontok[4], uu);

	vec2 b1[3];
	b1[0] = oszto_pont(b0[0], b0[1], uu);
	b1[1] = oszto_pont(b0[1], b0[2], uu);
	b1[2] = oszto_pont(b0[2], b0[3], uu);

	vec2 b2[2];
	b2[0] = oszto_pont(b1[0], b1[1], uu);
	b2[1] = oszto_pont(b1[1], b1[2], uu);

	vec2 b31 = oszto_pont(b2[0], b2[1], uu);

	return b31;
}

void de_Casteljau() {
	calculatePoints();

	vec2 b0[4];
	b0[0] = oszto_pont(points[6], points[7], u);
	b0[1] = oszto_pont(points[7], points[8], u);
	b0[2] = oszto_pont(points[8], points[9], u);
	b0[3] = oszto_pont(points[9], points[10], u);

	if (displayPoligon) {
		glLineWidth(1.0);
		glColor3f(0.0, 0.5, 0.5);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 4; i++) {
			glVertex2f(b0[i].x, b0[i].y);
		}
		glEnd();

		glPointSize(8.0);
		glColor3f(0.75, 0.5, 0.5);
		glBegin(GL_POINTS);
		for (int i = 0; i < 4; i++) {
			glVertex2f(b0[i].x, b0[i].y);
		}
		glEnd();
	}

	vec2 b1[3];
	b1[0] = oszto_pont(b0[0], b0[1], u);
	b1[1] = oszto_pont(b0[1], b0[2], u);
	b1[2] = oszto_pont(b0[2], b0[3], u);

	if (displayPoligon) {
		glLineWidth(1.0);
		glColor3f(0.0, 0.5, 0.5);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 3; i++) {
			glVertex2f(b1[i].x, b1[i].y);
		}
		glEnd();

		glPointSize(8.0);
		glColor3f(0.75, 0.5, 0.5);
		glBegin(GL_POINTS);
		for (int i = 0; i < 3; i++) {
			glVertex2f(b1[i].x, b1[i].y);
		}
		glEnd();
	}

	vec2 b2[2];
	b2[0] = oszto_pont(b1[0], b1[1], u);
	b2[1] = oszto_pont(b1[1], b1[2], u);

	if (displayPoligon) {
		glLineWidth(1.0);
		glColor3f(0.0, 0.5, 0.5);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 2; i++) {
			glVertex2f(b2[i].x, b2[i].y);
		}
		glEnd();

		glPointSize(8.0);
		glColor3f(0.75, 0.5, 0.5);
		glBegin(GL_POINTS);
		for (int i = 0; i < 2; i++) {
			glVertex2f(b2[i].x, b2[i].y);
		}
		glEnd();
	}


	vec2 b31 = oszto_pont(b2[0], b2[1], u);

	glPointSize(8);
	glColor3f(0.0, 0.0, 0.0);
	glBegin(GL_POINTS);
		glVertex2f(b31.x, b31.y);
	glEnd();

	std::vector<vec2> pp(5);
	pp[0] = points[6];
	pp[1] = points[7];
	pp[2] = points[8];
	pp[3] = points[9];
	pp[4] = points[10];

	glLineWidth(2);
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINE_STRIP);
	for (float t = 0; t <= 1; t += 0.001) {
		vec2 curvePoint = calculateCurvePoint(pp, t);
		glVertex2f(curvePoint.x, curvePoint.y);
	}
	glEnd();
}

void init() {
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glShadeModel(GL_SMOOTH);
	gluOrtho2D(0.0, winWidth, 0.0, winHeight);
}

void myDisplay() {
	glClear(GL_COLOR_BUFFER_BIT);

	hermite();
	bernstein();
	de_Casteljau();

	if (displayPoligon) {
		displayControlPolygon();
	}

	displayControlPoints();
	glutSwapBuffers();
}

GLint getActivePoint(std::vector<vec2> p, GLint size, GLint sens, GLint x, GLint y) {
	GLint i, s = sens * sens;
	vec2 P = { (float)x, (float)y };

	for (i = 0; i < size; i++)
		if (dist(p[i], P) < s)
			return i;
	return -1;
}

void processMouse(GLint button, GLint action, GLint xMouse, GLint yMouse) {
	GLint i;
	if (button == GLUT_LEFT_BUTTON && action == GLUT_DOWN)
		if ((i = getActivePoint(points, points.size(), pickRadius, xMouse, winHeight - yMouse)) != -1)
			dragged = i;
	if (button == GLUT_LEFT_BUTTON && action == GLUT_UP)
		dragged = -1;
}

void processMouseActiveMotion(GLint xMouse, GLint yMouse) {
	GLint i;
	if (dragged >= 0) {
		points[dragged].x = xMouse;
		points[dragged].y = winHeight - yMouse;
		if (points[dragged].x + pickRadius >= winWidth)
			points[dragged].x = winWidth-pickRadius;
		if (points[dragged].x - pickRadius <= 0)
			points[dragged].x = 0+pickRadius;
		if (points[dragged].y + pickRadius >= winHeight)
			points[dragged].y = winHeight-pickRadius;
		if (points[dragged].y - pickRadius <= 0)
			points[dragged].y = 0+pickRadius;
		glutPostRedisplay();
	}
}

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	//	esc key
	case 27:
		exit(0);
		break;
	case 'w':
		if (u < 1) u += 0.002;
		else u = 1;
		break;
	case 's':
		if (u > 0) u -= 0.002;
		else u = 0;
		break;
	case 'r':
		points = originalPoints;
		displayPoligon = true;
		break;
	case 'p':
			displayPoligon = !displayPoligon;
	}
	glutPostRedisplay();
}

int main(int argc, char** argv) {
	//	Init
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(winWidth, winHeight);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("1. Beadandó");

	//	Custom functions
	glutDisplayFunc(myDisplay);
	glutMouseFunc(processMouse);
	glutMotionFunc(processMouseActiveMotion);
	glutKeyboardFunc(keyboard);

	init();
	glutMainLoop();
	return 0;
}
