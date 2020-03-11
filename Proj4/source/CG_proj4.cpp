#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <deque>
#include <queue>
#include <GL/glut.h>

#define PRINT 0
#define PRINTSPECIAL 1
#define PRINT_VERTEX 0

#define MAX_WIN_SIZE 1000
#define MTM_SIZE 4
#define MAX_MODEL 20
#define MAX_MODEL_NUM 4000
#define MAX_POINT_IN_FACE 10
#define INF_FAR 999.0
#define MAX_LIGHT 10
using namespace std;

typedef float wvmtype;
int width;
int height;
float background_r = 0.0;
float background_g = 0.0;
float background_b = 0.0;
float alpha = 0.0;
float ambient_Kar = 0.0;
float ambient_Kag = 0.0;
float ambient_Kab = 0.0;

bool IsExit = false;

int num_ASCModel = 0;
int num_line = 0;
float Eye_Matirx[MTM_SIZE][MTM_SIZE];
float Project_Matrix[MTM_SIZE][MTM_SIZE];

float EyeX, EyeY, EyeZ;
float WzNear, WzFar, WhFOV;

wvmtype WVM[MTM_SIZE][MTM_SIZE];
struct frame_struct {
	float vxl;
	float vxr;
	float vyb;
	float vyt;
};
struct frame_struct frame;
struct Buffer_struct {
	float z;
	float r, g, b;
	Buffer_struct() {}
	Buffer_struct(float _z, float _r, float _g, float _b) {
		z = _z;
		r = _r;
		g = _g;
		b = _b;
	}
}Buffer[MAX_WIN_SIZE][MAX_WIN_SIZE];
float ModelingTransformMatrix[MTM_SIZE][MTM_SIZE]
= { { 1, 0, 0, 0 },
	{ 0, 1, 0, 0 },
	{ 0, 0, 1, 0 },
	{ 0, 0, 0, 1 } };
struct ASCModel_struct {
	int num_vertex;
	int num_face;
	float vertex[MAX_MODEL_NUM][MTM_SIZE];
	int face[MAX_MODEL_NUM][MAX_POINT_IN_FACE];
	float R, G, B;
	float Kd, Ks, N;
}ASCModel[MAX_MODEL];
struct PointWithColor {
	int y;
	float x, z;
	float r, g, b;
	PointWithColor() :x(-1), y(-1) {}
	friend bool operator<(const PointWithColor& p1, const PointWithColor& p2) {
		return (p1.x > p2.x);
	}
};
struct Light_struct {
	bool enable;
	float Ipr, Ipg, Ipb;
	float X, Y, Z;
	Light_struct() {
		enable = false;
	}
}Light[MAX_LIGHT];

ifstream fin;

void init(void)
{
	glClearColor(background_r, background_g, background_b, alpha);
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0, width, 0, height);
	glShadeModel(GL_FLAT);
	glFlush();
}

void initial() {
	Buffer_struct infBuffer(INF_FAR, background_r, background_g, background_b);
	Buffer_struct blackBuffer(INF_FAR, 0, 0, 0);
	cout << frame.vxl << frame.vxr << frame.vyb << frame.vyt << endl;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (j >= frame.vxl && j < frame.vxr && i >= frame.vyb && i < frame.vyt) {
				Buffer[i][j] = infBuffer;
			}
			else {
				Buffer[i][j] = blackBuffer;
			}
		}
	}
}
void setAmbient(float r, float g, float b) {
	ambient_Kar = r;
	ambient_Kag = g;
	ambient_Kab = b;
}
void setBackground(float _background_r, float _background_g, float _background_b) {
	background_r = _background_r;
	background_g = _background_g;
	background_b = _background_b;
}
void setLight(int _ID, float _Ipr, float _Ipg, float _Ipb, float _X, float _Y, float _Z) {
	Light[_ID].enable = true;
	Light[_ID].Ipr = _Ipr;
	Light[_ID].Ipg = _Ipg;
	Light[_ID].Ipb = _Ipb;
	Light[_ID].X = _X;
	Light[_ID].Y = _Y;
	Light[_ID].Z = _Z;
}

void Matrix_Multi_Matrix(float a[][MTM_SIZE], float b[][MTM_SIZE], float c[][MTM_SIZE]) {
	float temp[MTM_SIZE][MTM_SIZE];

	float multi_sum = 0;
	for (int i = 0; i < MTM_SIZE; i++) {
		for (int j = 0; j < MTM_SIZE; j++) {
			multi_sum = 0;
			for (int k = 0; k < MTM_SIZE; k++) {
				multi_sum += b[k][j] * a[i][k];
			}
			temp[i][j] = multi_sum;
		}
	}

	// Fix to avoid input of Matrix_A is input of Matrix C
	for (int i = 0; i < MTM_SIZE; i++) {
		for (int j = 0; j < MTM_SIZE; j++) {
			c[i][j] = temp[i][j];
		}
	}
}
template<class T>
void printMatrix(T* matrix, int row, int col, string str) {
	cout << endl;
	cout << "Matrix <<< " << str << " >>> = " << endl;
	for (int i = 0; i < row; i++) {
		cout << "\t[";
		for (int j = 0; j < col; j++) {
			cout.setf(ios::fixed);
			cout << setw(5) << setprecision(2) << matrix[i][j] << ", ";
		}
		cout << " ]" << endl;
	}
}
template<class T>
void update_Matrix(T resource[][MTM_SIZE], T destination[][MTM_SIZE]) {
	for (int i = 0; i < MTM_SIZE; i++) {
		for (int j = 0; j < MTM_SIZE; j++) {
			destination[i][j] = resource[i][j];
		}
	}
}
void UnitizeVector(float vector[]) {
	float diver = sqrt(vector[0] * vector[0] +
		vector[1] * vector[1] +
		vector[2] * vector[2]);
	for (int i = 0; i < 3; i++) {
		vector[i] /= diver;
	}
}
void Cross_Multi(float destination[], float a[], float b[]) {
	destination[0] = a[1] * b[2] - a[2] * b[1];
	destination[1] = a[2] * b[0] - a[0] * b[2];
	destination[2] = a[0] * b[1] - a[1] * b[0];
}
void getNormalVector(float* destination, int objectIndex, int faceIndex) {
	float point1[3];
	float point2[3];
	float point3[3];
	float vector1[3];
	float vector2[3];
	for (int i = 0; i < 3; i++) {
		point1[i] = ASCModel[objectIndex].vertex[ASCModel[objectIndex].face[faceIndex][1] - 1][i];
		point2[i] = ASCModel[objectIndex].vertex[ASCModel[objectIndex].face[faceIndex][2] - 1][i];
		point3[i] = ASCModel[objectIndex].vertex[ASCModel[objectIndex].face[faceIndex][3] - 1][i];
		vector1[i] = point2[i] - point1[i];
		vector2[i] = point3[i] - point1[i];
	}
	Cross_Multi(destination, vector1, vector2);
	for (int i = 0; i < 3; i++) {
		destination[i] *= -1;
	}
}

float ReluDot_Multi(float a[], float b[]) {
	float sum = 0;
	for (int i = 0; i < 3; i++) {
		sum += a[i] * b[i];
	}
	//return sum;
	return (sum > 0) ? sum : 0; // Maybe Wrong
}

float interpolation(float p1, float p2, float q1, float q2, float q) {
	return ((q - q1) / (q2 - q1) * (p2 - p1) + p1);
}

float getColor(int objectIndex, int faceIndex, int pointIndex, char whichColor) {
	//cout << objectIndex << " " << faceIndex << " " << pointIndex << endl;
	float sumColor, ambientr, ambientg, ambientb, diffuse_r, diffuse_g, diffuse_b, specular_r, specular_g, specular_b;
	ambientr = ambient_Kar;
	ambientg = ambient_Kag;
	ambientb = ambient_Kab;

	diffuse_r = 0;
	diffuse_g = 0;
	diffuse_b = 0;

	specular_r = 0;
	specular_g = 0;
	specular_b = 0;

	float normalVector[3];
	getNormalVector(normalVector, objectIndex, faceIndex);
	UnitizeVector(normalVector);

	for (int i = 0; i < MAX_LIGHT; i++) {
		if (Light[i].enable) {
			float lightVector[3];
			lightVector[0] = Light[i].X - ASCModel[objectIndex].vertex[pointIndex][0];
			lightVector[1] = Light[i].Y - ASCModel[objectIndex].vertex[pointIndex][1];
			lightVector[2] = Light[i].Z - ASCModel[objectIndex].vertex[pointIndex][2];
			UnitizeVector(lightVector);
			float NdotL = ReluDot_Multi(normalVector, lightVector);
			float reflectVector[3];
			reflectVector[0] = 2 * NdotL * normalVector[0] - lightVector[0];
			reflectVector[1] = 2 * NdotL * normalVector[1] - lightVector[1];
			reflectVector[2] = 2 * NdotL * normalVector[2] - lightVector[2];
			UnitizeVector(reflectVector);
			float viewVector[3];
			viewVector[0] = EyeX - ASCModel[objectIndex].vertex[pointIndex][0];
			viewVector[1] = EyeY - ASCModel[objectIndex].vertex[pointIndex][1];
			viewVector[2] = EyeZ - ASCModel[objectIndex].vertex[pointIndex][2];
			UnitizeVector(viewVector);
			diffuse_r += ASCModel[objectIndex].Kd * ReluDot_Multi(normalVector, lightVector)* Light[i].Ipr;
			diffuse_g += ASCModel[objectIndex].Kd * ReluDot_Multi(normalVector, lightVector)* Light[i].Ipg;
			diffuse_b += ASCModel[objectIndex].Kd * ReluDot_Multi(normalVector, lightVector)* Light[i].Ipb;

			float specular_last = pow(ReluDot_Multi(reflectVector, viewVector), ASCModel[objectIndex].N);
			specular_r += ASCModel[objectIndex].Ks * specular_last * Light[i].Ipr;
			specular_g += ASCModel[objectIndex].Ks * specular_last * Light[i].Ipg;
			specular_b += ASCModel[objectIndex].Ks * specular_last * Light[i].Ipb;
		}
	}

	if (whichColor == 'r')
		sumColor = (ambientr + diffuse_r) * ASCModel[objectIndex].R + specular_r;
	else if (whichColor == 'g')
		sumColor = (ambientg + diffuse_g) * ASCModel[objectIndex].G + specular_g;
	else
		sumColor = (ambientb + diffuse_b) * ASCModel[objectIndex].B + specular_b;

	return sumColor;
}
void reset() {
	for (int i = 0; i < MTM_SIZE; i++) {
		for (int j = 0; j < MTM_SIZE; j++) {
			if (i == j)
				ModelingTransformMatrix[i][j] = 1.0;
			else
				ModelingTransformMatrix[i][j] = 0;
		}
	}
#if PRINT
	printMatrix(ModelingTransformMatrix, 4, 4, "ModelingTransformMatrix");
#endif
}
void clearData() {
	num_line = 0;
}
void clearScreen() {
	glClearColor(background_r, background_g, background_b, alpha);
	glClear(GL_COLOR_BUFFER_BIT);
	glFlush();
}

void Eye_Transform(float PX, float PY, float PZ, float CX, float CY, float CZ, float Tilt) {
	float Eye_Translation_Matrix[MTM_SIZE][MTM_SIZE]
		= { { 1, 0, 0, -PX },
			{ 0, 1, 0, -PY },
			{ 0, 0, 1, -PZ },
			{ 0, 0, 0, 1 }
	};
#if PRINT
	printMatrix(Eye_Translation_Matrix, 4, 4, "Eye_Translation_Matrix");
#endif
	float Eye_Mirror_Matrix[MTM_SIZE][MTM_SIZE]
		= { { -1, 0, 0, 0 },
			{ 0, 1, 0, 0 },
			{ 0, 0, 1, 0 },
			{ 0, 0, 0, 1 }
	};
#if PRINT
	printMatrix(Eye_Mirror_Matrix, 4, 4, "Eye_Mirror_Matrix");
#endif

	float TiltDegree = Tilt * M_PI / 180.0;
	float Eye_Tilt_Matrix[MTM_SIZE][MTM_SIZE]
		= { { cos(TiltDegree), sin(TiltDegree), 0, 0 },
			{ -sin(TiltDegree), cos(TiltDegree), 0, 0 },
			{ 0, 0, 1, 0 },
			{ 0, 0, 0, 1 } };
#if PRINT
	printMatrix(Eye_Tilt_Matrix, 4, 4, "Eye_Tilt_Matrix");
#endif

	float view_vector[3] = { CX - PX, CY - PY, CZ - PZ };
	float top_vector[3] = { 0, 1, 0 };

	float vector3[3];
	memcpy(vector3, view_vector, sizeof(view_vector));
	UnitizeVector(vector3);
	float vector1[3];
	Cross_Multi(vector1, top_vector, view_vector);
	UnitizeVector(vector1);
	float vector2[3];
	Cross_Multi(vector2, vector3, vector1);
	UnitizeVector(vector2);

	float GRM[MTM_SIZE][MTM_SIZE]
		= { { vector1[0], vector1[1], vector1[2], 0 },
			{ vector2[0], vector2[1], vector2[2], 0 },
			{ vector3[0], vector3[1], vector3[2], 0 },
			{ 0, 0, 0, 1 } };
#if PRINT
	printMatrix(GRM, 4, 4, "GRM");
#endif

	float new_Eye_Matrix[MTM_SIZE][MTM_SIZE]
		= { { 1, 0, 0, 0 },
			{ 0, 1, 0, 0 },
			{ 0, 0, 1, 0 },
			{ 0, 0, 0, 1 } };
	Matrix_Multi_Matrix(Eye_Translation_Matrix, new_Eye_Matrix, new_Eye_Matrix);
	Matrix_Multi_Matrix(GRM, new_Eye_Matrix, new_Eye_Matrix);
	Matrix_Multi_Matrix(Eye_Mirror_Matrix, new_Eye_Matrix, new_Eye_Matrix);
	Matrix_Multi_Matrix(Eye_Tilt_Matrix, new_Eye_Matrix, new_Eye_Matrix);
	update_Matrix(new_Eye_Matrix, Eye_Matirx);
#if PRINT
	printMatrix(Eye_Matirx, 4, 4, "Eye_Matrix");
#endif
}
void Project_Transform(float zNear, float zFar, float hFOV) {
	float PM4_3 = tan(hFOV * M_PI / 180.0);
	float PM3_3 = (zFar / (zFar - zNear)) * PM4_3;
	float PM3_4 = ((zNear * zFar) / (zNear - zFar)) * PM4_3;
	float new_Project_Matrix[MTM_SIZE][MTM_SIZE]
		= { { 1, 0, 0, 0 },
			{ 0, 1, 0, 0 },
			{ 0, 0, PM3_3, PM3_4 },
			{ 0, 0, PM4_3, 0 } };
	update_Matrix(new_Project_Matrix, Project_Matrix);
}

void drawColorDot(int x, int y, float c_r, float c_g, float c_b) {
	if (x <= frame.vxl || x > frame.vxr || y <= frame.vyb || y > frame.vyt) {
		c_r = 0.0;
		c_g = 0.0;
		c_b = 0.0;
	}
	glBegin(GL_POINTS);
	// set the color of dot
	glColor3f(c_r, c_g, c_b);
	glVertex2i(x, y);

	glEnd();
}
void redraw() {
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			drawColorDot(j, i, Buffer[i][j].r, Buffer[i][j].g, Buffer[i][j].b);
		}
	}
	glFlush();
}

void observer(float PX, float PY, float PZ, float CX, float CY, float CZ,
	float Tilt, float zNear, float zFar, float hFOV) {
	Eye_Transform(PX, PY, PZ, CX, CY, CZ, Tilt);

	Project_Transform(zNear, zFar, hFOV);

	EyeX = PX;
	EyeY = PY;
	EyeZ = PZ;

	WzNear = zNear;
	WzFar = zFar;
	WhFOV = hFOV;
}

void object(string objectname, float _R, float _G, float _B, float _Kd, float _Ks, float _N) {
	ifstream fin(objectname);
	if (fin.is_open()) {
		cout << "\topen the " << objectname << " successfully" << endl;
	}
	else {
		cout << "\tCan't open the " << objectname << endl;
		return;
	}
	// get number of vertex and face first
	fin >> ASCModel[num_ASCModel].num_vertex;
	fin >> ASCModel[num_ASCModel].num_face;
	// read vertex one by one
	for (int i = 0; i < ASCModel[num_ASCModel].num_vertex; i++) {
		for (int j = 0; j < 3; j++) {
			fin >> ASCModel[num_ASCModel].vertex[i][j];
		}
		ASCModel[num_ASCModel].vertex[i][3] = 1;//fix to multiplication
	}

	// read face one by one
	for (int i = 0; i < ASCModel[num_ASCModel].num_face; i++) {
		fin >> ASCModel[num_ASCModel].face[i][0];
		for (int j = 1; j <= ASCModel[num_ASCModel].face[i][0]; j++) {
			fin >> ASCModel[num_ASCModel].face[i][j];
		}
	}

	// Convert Object-Space to World-Space
	float multi_sum = 0;
	float WorldSpace[MTM_SIZE];
	for (int i = 0; i < ASCModel[num_ASCModel].num_vertex; i++) {
		for (int j = 0; j < MTM_SIZE; j++) {
			multi_sum = 0;
			for (int k = 0; k < MTM_SIZE; k++) {
				multi_sum += ASCModel[num_ASCModel].vertex[i][k]
					* ModelingTransformMatrix[j][k];
			}
			WorldSpace[j] = multi_sum;
		}
		// Update(cover) the Object-Space
		for (int j = 0; j < MTM_SIZE; j++) {
			ASCModel[num_ASCModel].vertex[i][j] = WorldSpace[j];
		}
	}

	ASCModel[num_ASCModel].R = _R;
	ASCModel[num_ASCModel].G = _G;
	ASCModel[num_ASCModel].B = _B;
	ASCModel[num_ASCModel].Kd = _Kd;
	ASCModel[num_ASCModel].Ks = _Ks;
	ASCModel[num_ASCModel].N = _N;

	num_ASCModel++;
}
void viewport(float vxl, float vxr, float vyb, float vyt) {
	float Vlength = vxr - vxl;
	float Vheight = vyt - vyb;
	Project_Matrix[1][1] = Vlength / Vheight;// Set Ar
#if PRINT
	printMatrix(Project_Matrix, 4, 4, "Project_Matrix");
#endif

	float wxl, wxr, wyb, wyt;
	wxl = wyb = -1.0;
	wxr = wyt = 1.0;
	wvmtype scaling_x = (wvmtype)((vxr - vxl) / (wxr - wxl));
	wvmtype scaling_y = (wvmtype)((vyt - vyb) / (wyt - wyb));
	wvmtype shift_x = (wvmtype)(vxl - scaling_x * wxl);
	wvmtype shift_y = (wvmtype)(vyb - scaling_y * wyb);
	scaling_x = (scaling_x)* width / 2;
	scaling_y = (scaling_y)* height / 2;
	shift_x = (shift_x + 1) * width / 2;
	shift_y = (shift_y + 1) * height / 2;

	frame.vxl = (int)((vxl + 1) * width / 2);
	frame.vxr = (int)((vxr + 1) * width / 2);
	frame.vyb = (int)((vyb + 1) * height / 2);
	frame.vyt = (int)((vyt + 1) * height / 2);
	wvmtype new_WVM_Matrix[MTM_SIZE][MTM_SIZE]
		= { { scaling_x, 0, shift_x, 0 },
			{ 0, scaling_y, shift_y, 0 },
			{ 0, 0, 1, 0 },
			{ 0, 0, 0, 1 } };
	update_Matrix(new_WVM_Matrix, WVM);
#if PRINTSPECIAL
	printMatrix(WVM, 4, 4, "WVM Matrix");
#endif
}

void scale(float sx, float sy, float sz) {
	float scaling_matrix[MTM_SIZE][MTM_SIZE]
		= { { sx, 0, 0, 0 },
			{ 0, sy, 0, 0 },
			{ 0, 0, sz, 0 },
			{ 0, 0, 0, 1 }
	};

	float new_matrix[MTM_SIZE][MTM_SIZE];

	Matrix_Multi_Matrix(scaling_matrix, ModelingTransformMatrix, new_matrix);

	update_Matrix(new_matrix, ModelingTransformMatrix);
#if PRINT
	printMatrix(ModelingTransformMatrix, 4, 4, "ModelingTransformMatrix");
#endif
}
void translate(float tx, float ty, float tz) {
	float translate_matrix[MTM_SIZE][MTM_SIZE]
		= { { 1, 0, 0, tx },
			{ 0, 1, 0, ty },
			{ 0, 0, 1, tz },
			{ 0, 0, 0, 1 }
	};

	float new_matrix[MTM_SIZE][MTM_SIZE];

	Matrix_Multi_Matrix(translate_matrix, ModelingTransformMatrix, new_matrix);

	update_Matrix(new_matrix, ModelingTransformMatrix);
#if PRINT
	printMatrix(ModelingTransformMatrix, 4, 4, "ModelingTransformMatrix");
#endif
}
void rotate(float Xdegree, float Ydegree, float Zdegree) {
	float new_matrix[MTM_SIZE][MTM_SIZE];

	if (Xdegree != 0) {
		float rad = Xdegree * M_PI / 180.0;
		float rotation_matrix[MTM_SIZE][MTM_SIZE]
			= { { 1, 0, 0, 0 },
				{ 0, cos(rad), -sin(rad), 0 },
				{ 0, sin(rad), cos(rad), 0 },
				{ 0, 0, 0, 1 } };

		Matrix_Multi_Matrix(rotation_matrix, ModelingTransformMatrix, new_matrix);
		update_Matrix(new_matrix, ModelingTransformMatrix);
#if PRINT
		printMatrix(ModelingTransformMatrix, 4, 4, "ModelingTransformMatrix");
#endif
	}

	if (Ydegree != 0) {
		float rad = Ydegree * M_PI / 180.0;
		float rotation_matrix[MTM_SIZE][MTM_SIZE]
			= { { cos(rad), 0, sin(rad), 0 },
				{ 0, 1, 0, 0 },
				{ -sin(rad), 0, cos(rad), 0 },
				{ 0, 0, 0, 1 } };

		Matrix_Multi_Matrix(rotation_matrix, ModelingTransformMatrix, new_matrix);
		update_Matrix(new_matrix, ModelingTransformMatrix);
#if PRINT
		printMatrix(ModelingTransformMatrix, 4, 4, "ModelingTransformMatrix");
#endif
	}

	if (Zdegree != 0) {
		float rad = Zdegree * M_PI / 180.0;
		float rotation_matrix[MTM_SIZE][MTM_SIZE]
			= { { cos(rad), -sin(rad), 0, 0 },
				{ sin(rad), cos(rad), 0, 0 },
				{ 0, 0, 1, 0 },
				{ 0, 0, 0, 1 } };

		Matrix_Multi_Matrix(rotation_matrix, ModelingTransformMatrix, new_matrix);
		update_Matrix(new_matrix, ModelingTransformMatrix);
#if PRINT
		printMatrix(ModelingTransformMatrix, 4, 4, "ModelingTransformMatrix");
#endif
	}
}

void drawFace(const deque<PointWithColor>& colorPointVector, int min_y, int max_y) {
	deque< priority_queue<PointWithColor> > vecPriQue;
	for (int i = min_y; i <= max_y; i++) {
		priority_queue<PointWithColor> empty;
		vecPriQue.push_back(empty);
	}
	for (int i = 0; i < colorPointVector.size(); i++) {
		int i_next = (i + 1 == colorPointVector.size()) ? 0 : i + 1;
		//int a = abs(colorPointVector[i].x - colorPointVector[i_next].x);
		int b = abs(colorPointVector[i].y - colorPointVector[i_next].y);
		//int max_dis = a > b ? a : b;
		int max_dis = b;
		//cout << a << "  " << b << "  " << max_dis << endl;
		if (max_dis != 0) {
			PointWithColor newColorPoint;
			for (int j = 0; j < max_dis; j++) {
				newColorPoint.x = colorPointVector[i].x + j
					* (colorPointVector[i_next].x - colorPointVector[i].x) / max_dis;
				newColorPoint.y = colorPointVector[i].y + j
					* (colorPointVector[i_next].y - colorPointVector[i].y) / max_dis;
				//cout << newColorPoint.x << "   " << newColorPoint.y << endl;
				newColorPoint.z = interpolation(colorPointVector[i].z, colorPointVector[i_next].z,
					colorPointVector[i].y, colorPointVector[i_next].y, newColorPoint.y);
				//cout << newColorPoint.z << endl;
				newColorPoint.r = interpolation(colorPointVector[i].r, colorPointVector[i_next].r,
					colorPointVector[i].y, colorPointVector[i_next].y, newColorPoint.y);
				newColorPoint.g = interpolation(colorPointVector[i].g, colorPointVector[i_next].g,
					colorPointVector[i].y, colorPointVector[i_next].y, newColorPoint.y);
				newColorPoint.b = interpolation(colorPointVector[i].b, colorPointVector[i_next].b,
					colorPointVector[i].y, colorPointVector[i_next].y, newColorPoint.y);
				vecPriQue[newColorPoint.y - min_y].push(newColorPoint);
				//drawColorDot(newColorPoint.x, newColorPoint.y, newColorPoint.r, newColorPoint.g, newColorPoint.b);
			}
		}
		else {
			vecPriQue[colorPointVector[i].y - min_y].push(colorPointVector[i]);
			vecPriQue[colorPointVector[i_next].y - min_y].push(colorPointVector[i_next]);
		}
	}
	for (int i = 0; i < vecPriQue.size(); i++) {
		if (!vecPriQue[i].empty()) {
			PointWithColor left = vecPriQue[i].top();
			while (vecPriQue[i].size() > 1)
				vecPriQue[i].pop();
			PointWithColor right = vecPriQue[i].top();
			vecPriQue[i].pop();
			if (right.x - left.x == 0) {
				if (left.z < Buffer[left.y][(int)left.x].z) {
					Buffer[left.y][(int)left.x].z = left.z;
					Buffer[left.y][(int)left.x].r = left.r;
					Buffer[left.y][(int)left.x].g = left.g;
					Buffer[left.y][(int)left.x].b = left.b;
				}
				continue;
			}

			int step_x = left.x - 1;
			float step_z = (right.z - left.z) / (right.x - left.x);     float depth_z = left.z - step_z;
			float step_r = (right.r - left.r) / (right.x - left.x);     float cur_r = left.r - step_r;
			float step_g = (right.g - left.g) / (right.x - left.x);     float cur_g = left.g - step_g;
			float step_b = (right.b - left.b) / (right.x - left.x);     float cur_b = left.b - step_b;

			int steps = right.x - left.x;
			for (int j = 0; j <= steps; j++) {
				step_x++;
				depth_z += step_z;
				cur_r += step_r;    cur_g += step_g;    cur_b += step_b;
				if (depth_z < Buffer[left.y][step_x].z) {
					Buffer[left.y][step_x].z = depth_z;
					Buffer[left.y][step_x].r = cur_r;
					Buffer[left.y][step_x].g = cur_g;
					Buffer[left.y][step_x].b = cur_b;
				}
			}
		}
	}
}

void display() {
	//viewport(0, width, 0, height);

	initial();
	//clearScreen();
	//clearData();

	deque<PointWithColor> colorPointVector;
	PointWithColor newColorPoint;

	float multi_sum = 0;
	float Matrix[MTM_SIZE];
	ASCModel_struct Cur_ASCModel;
	for (int index = 0; index < num_ASCModel; index++) {
		Cur_ASCModel = ASCModel[index];
#if PRINT
		cout << "World-Space:" << endl;
		cout << Cur_ASCModel.vertex[PRINT_VERTEX][0] << " "
			<< Cur_ASCModel.vertex[PRINT_VERTEX][1] << " "
			<< Cur_ASCModel.vertex[PRINT_VERTEX][2] << " "
			<< Cur_ASCModel.vertex[PRINT_VERTEX][3] << " " << endl;
#endif
		// Convert World-Space to Eye-Space
		for (int i = 0; i < ASCModel[index].num_vertex; i++) {
			for (int j = 0; j < MTM_SIZE; j++) {
				multi_sum = 0;
				for (int k = 0; k < MTM_SIZE; k++) {
					multi_sum += Cur_ASCModel.vertex[i][k]
						* Eye_Matirx[j][k];
				}
				Matrix[j] = multi_sum;
			}
			// Update(cover) the Object-Space
			for (int j = 0; j < MTM_SIZE; j++) {
				Cur_ASCModel.vertex[i][j] = Matrix[j];
			}
		}
#if PRINT
		cout << "Eye-Space:" << endl;
		cout << Cur_ASCModel.vertex[PRINT_VERTEX][0] << " "
			<< Cur_ASCModel.vertex[PRINT_VERTEX][1] << " "
			<< Cur_ASCModel.vertex[PRINT_VERTEX][2] << " "
			<< Cur_ASCModel.vertex[PRINT_VERTEX][3] << " " << endl;
#endif
		// Convert Eye-Space to Project-Space
		for (int i = 0; i < ASCModel[index].num_vertex; i++) {
			for (int j = 0; j < MTM_SIZE; j++) {
				multi_sum = 0;
				for (int k = 0; k < MTM_SIZE; k++) {
					multi_sum += Cur_ASCModel.vertex[i][k]
						* Project_Matrix[j][k];
				}
				Matrix[j] = multi_sum;
			}
			// Update(cover) the Object-Space
			for (int j = 0; j < MTM_SIZE; j++) {
				Cur_ASCModel.vertex[i][j] = Matrix[j] / Matrix[3];
			}
			//Cur_ASCModel.vertex[i][2] = 1;//Fix to three-dimensional homogeneous
		}
#if PRINT
		cout << "Project-Space:" << endl;
		cout << Cur_ASCModel.vertex[PRINT_VERTEX][0] << " "
			<< Cur_ASCModel.vertex[PRINT_VERTEX][1] << " "
			<< Cur_ASCModel.vertex[PRINT_VERTEX][2] << " "
			<< Cur_ASCModel.vertex[PRINT_VERTEX][3] << " " << endl;
#endif
		// Convert Project-Space to Screen-Space
		for (int i = 0; i < ASCModel[index].num_vertex; i++) {
			for (int j = 0; j < MTM_SIZE; j++) {
				multi_sum = 0;
				for (int k = 0; k < MTM_SIZE; k++) {
					multi_sum += Cur_ASCModel.vertex[i][k]
						* WVM[j][k];
				}
				Matrix[j] = multi_sum;
			}
			// Update(cover) the Object-Space
			for (int j = 0; j < MTM_SIZE; j++) {
				Cur_ASCModel.vertex[i][j] = Matrix[j];
			}
		}
#if PRINT
		cout << "Screen-Space:" << endl;
		cout << Cur_ASCModel.vertex[PRINT_VERTEX][0] << " "
			<< Cur_ASCModel.vertex[PRINT_VERTEX][1] << " "
			<< Cur_ASCModel.vertex[PRINT_VERTEX][2] << " "
			<< Cur_ASCModel.vertex[PRINT_VERTEX][3] << " " << endl;
#endif
		// Display
		for (int i = 0; i < Cur_ASCModel.num_face; i++) {
			colorPointVector.clear();
			int min_y = INT_MAX;
			int max_y = INT_MIN;
			for (int j = 1; j <= Cur_ASCModel.face[i][0]; j++) {
				newColorPoint.x = Cur_ASCModel.vertex[Cur_ASCModel.face[i][j] - 1][0];
				newColorPoint.y = Cur_ASCModel.vertex[Cur_ASCModel.face[i][j] - 1][1];
				if (newColorPoint.y < min_y)
					min_y = newColorPoint.y;
				if (newColorPoint.y > max_y)
					max_y = newColorPoint.y;
				newColorPoint.z = Cur_ASCModel.vertex[Cur_ASCModel.face[i][j] - 1][2];
				newColorPoint.r = getColor(index, i, Cur_ASCModel.face[i][j] - 1, 'r');
				newColorPoint.g = getColor(index, i, Cur_ASCModel.face[i][j] - 1, 'g');
				newColorPoint.b = getColor(index, i, Cur_ASCModel.face[i][j] - 1, 'b');
				colorPointVector.push_back(newColorPoint);
			}
			//sortColorPointVector(colorPointVector);
			drawFace(colorPointVector, min_y, max_y);
		}
		//num_line++;
	}

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			drawColorDot(j, i, Buffer[i][j].r, Buffer[i][j].g, Buffer[i][j].b);
		}
	}
	glFlush();
}
void readFile(bool& IsExit) {
	int _ID;
	float sx, sy, sz;
	float Xdegree, Ydegree, Zdegree;
	float tx, ty, tz;
	float vxl, vxr, vyt, vyb;
	float PX, PY, PZ, CX, CY, CZ, Tilt, zNear, zFar, hFOV;
	float AmKar, AmKag, AmKab, b_r, b_g, b_b;
	float _R, _G, _B, _Kd, _Ks, _N;
	float _Ipr, _Ipg, _Ipb, _X, _Y, _Z;
	string command, comment, objectname;

	while (!fin.eof()) {
		fin >> command;
		cout << "COMMAND -- " << command << " : " << endl;
		if (command == "scale") {
			fin >> sx >> sy >> sz;
			scale(sx, sy, sz);
		}
		else if (command == "rotate") {
			fin >> Xdegree >> Ydegree >> Zdegree;
			rotate(Xdegree, Ydegree, Zdegree);
		}
		else if (command == "translate") {
			fin >> tx >> ty >> tz;
			translate(tx, ty, tz);
		}
		else if (command == "reset") {
			reset();
		}
		else if (command == "object") {
			fin >> objectname;
			fin >> _R >> _G >> _B;
			fin >> _Kd >> _Ks >> _N;
			object(objectname, _R, _G, _B, _Kd, _Ks, _N);
		}
		else if (command == "observer") {
			fin >> PX >> PY >> PZ;
			fin >> CX >> CY >> CZ;
			fin >> Tilt;
			fin >> zNear >> zFar >> hFOV;
			observer(PX, PY, PZ, CX, CY, CZ, Tilt, zNear, zFar, hFOV);
		}
		else if (command == "viewport") {
			fin >> vxl >> vxr >> vyb >> vyt;
			viewport(vxl, vxr, vyb, vyt);
		}
		else if (command == "display") {
			display();
			cout << "\tObjects display successfully" << endl;
			system("pause");
		}
		else if (command == "clearData") {
			clearData();
			cout << "Data is cleared" << endl;
		}
		else if (command == "clearScreen") {
			clearScreen();
			cout << "Screen is cleared" << endl;
		}
		else if (command == "ambient") {
			fin >> AmKar >> AmKag >> AmKab;
			setAmbient(AmKar, AmKag, AmKab);
		}
		else if (command == "background") {
			fin >> b_r >> b_g >> b_b;
			setBackground(b_r, b_g, b_b);
		}
		else if (command == "light") {
			fin >> _ID;
			fin >> _Ipr >> _Ipg >> _Ipb >> _X >> _Y >> _Z;
			setLight(_ID - 1, _Ipr, _Ipg, _Ipb, _X, _Y, _Z);
		}
		else if (command == "end") {
			IsExit = true;
			fin.close();
			exit(0);
			return;
		}
		else if (command == "#") {
			getline(fin, comment);
			cout << "\t" << comment << endl;
		}
		cout << endl;
	}
}
void displayFunc(void) {
	// clear the entire window to the background color
	glClear(GL_COLOR_BUFFER_BIT);
	while (!IsExit) {
		glClearColor(background_r, background_g, background_b, alpha);
		readFile(IsExit);
		return;
	}

	if (IsExit) {
		redraw();
	}
	return;
}

int main(int argc, char **argv) {
	// init GLUT and create Window
	string filename = argc >= 2 ? argv[1] : "Lab4D.in";
	fin = ifstream(filename);
	system("pause");
	if (fin.is_open()) {
		cout << "open the file successfully" << endl;
	}
	else {
		cout << "Can't not open the file" << endl;
		return 0;
	}
	fin >> width >> height;
	cout << "width :" << width << " height: " << height << endl;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
	glutInitWindowSize(width, height);
	glutInitWindowPosition(100, 10);
	glutCreateWindow("Assignment4");
	glutDisplayFunc(displayFunc);
	init();
	glutMainLoop();
	return 0;
}