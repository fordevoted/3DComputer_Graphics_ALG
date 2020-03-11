#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <GL/glut.h>
#include <iomanip>
#include <deque>

using namespace std;

#define MTM_SIZE 4
#define MAX_MODEL 15
#define MAX_MODEL_NUM 4000

#define WMinusX 0
#define WPlusX 1
#define WMinusY 2
#define WPlusY 3
#define WMinusZ 4
#define Z0 5

#define R 0.2512
#define G 0.8455
#define B 0.5416

struct ASCModel_struct {
	int num_vertex;
	int num_face;
	float vertex[MAX_MODEL_NUM][MTM_SIZE];
	int face[MAX_MODEL_NUM][MAX_MODEL];
}ASCModel[MAX_MODEL];
typedef struct ASCModel_struct ASCModel_struct;
typedef float wvmtype;
struct Point {
	float x, y;
	float z;
	float w;
	int vertexNumber;
};
struct response_type {
	deque<Point> polygon;
	bool isClip;
};
typedef struct response_type responseType;

struct frame_struct {
	float vxl;
	float vxr;
	float vyb;
	float vyt;
};

struct frame_struct frame;
int width;
int height;
int num_line = 0;
int num_ASCModel = 0;
float WzNear, WzFar, WhFOV;
float view_vactor[3];
bool NoBackFace = false;

ASCModel_struct Cur_ASCModel;
ifstream fin;

float Eye_Matirx[MTM_SIZE][MTM_SIZE];
float Project_Matrix[MTM_SIZE][MTM_SIZE];
wvmtype WVM[MTM_SIZE][MTM_SIZE];
float ModelingTransformMatrix[MTM_SIZE][MTM_SIZE]
= { { 1, 0, 0, 0 },
	{ 0, 1, 0, 0 },
	{ 0, 0, 1, 0 },
	{ 0, 0, 0, 1 } };

int sign(int a) {
	return a > 0 ? 1 : (a < 0 ? -1 : 0);
}
void drawLine(int x1, int y1, int x2, int y2, float r, float g, float b) {
	int i, ds, dx, dy, x, y;
	float slope;
	glColor3f(r, g, b);
	glBegin(GL_POINTS);
	if (y1 == y2)         /* for a horizontal line */
	{
		x = min(x1, x2);
		y = y1;
		while (x < max(x1, x2))
		{
			glVertex2f(x, y);
			x += 1;
		}
	}
	else if (x1 == x2)  /* for a vertical line */
	{
		y = min(y1, y2);
		x = x1;
		while (y < max(y1, y2))
		{
			glVertex2f(x, y);
			y += 1;
		}
	}
	else if (abs(x2 - x1) == abs(y2 - y1))
	{
		x = x1;
		y = y1;
		do
		{
			glVertex2f(x, y);
			x += sign(x2 - x1);
			y += sign(y2 - y1);
		} while (x != x2 && y != y2);
		glVertex2f(x, y);
	}
	else                 /* for straight lines with other slopes, best solution with slope +/-1 */
	{
		slope = (float)(y2 - y1) / (float)(x2 - x1);
		dx = abs(x2 - x1);
		dy = abs(y2 - y1);
		ds = 2 * dy - dx;
		x = min(x1, x2);
		y = min(y1, y2);
		if (slope > 0 && slope < 1)
		{
			glVertex2f(x, y);
			while (x < max(x1, x2))
			{
				if (ds <= 0)
					ds += 2 * dy;
				else
				{
					ds += 2 * (dy - dx);
					y += 1;
				}
				x += 1;
				glVertex2f(x, y);
			}
		}
		else if (slope > 1)
		{
			glVertex2f(x, y);
			while (y < max(y1, y2))
			{
				if (ds <= 0)
					ds += 2 * dx;
				else
				{
					ds += 2 * (dx - dy);
					x += 1;
				}
				y += 1;
				glVertex2f(x, y);
			}
		}
		else if (slope<0 && slope>-1)
		{
			ds = 2 * dy - dx;
			y = max(y1, y2);
			glVertex2f(x, y);
			while (x < max(x1, x2))
			{
				if (ds <= 0)
					ds += 2 * dy;
				else
				{
					ds += 2 * (dy - dx);
					y -= 1;
				}
				x += 1;
				glVertex2f(x, y);
			}
		}
		else
		{
			ds = 2 * dx - dy;
			x = max(x1, x2);
			y = min(y1, y2);
			glVertex2f(x, y);
			while (y < max(y1, y2))
			{
				if (ds <= 0)
					ds += 2 * dx;
				else
				{
					ds += 2 * (dx - dy);
					x -= 1;
				}
				y += 1;
				glVertex2f(x, y);
			}
		}
	}
	glEnd();
}
void init(void)
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0, width, 0, height);

	/* initialize viewing values */
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glFlush();
}
template<class T>
void printMatrix(T* matrix, int row, int col, string message) {
	cout << endl;
	cout << "Matrix <<< " << message << " >>> = " << endl;
	for (int i = 0; i < row; i++) {
		cout << "\t[";
		for (int j = 0; j < col; j++) {
			cout.setf(ios::fixed);
			cout << setw(5) << setprecision(2) << matrix[i][j] << ", ";
		}
		cout << " ]" << endl;
	}
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
	printMatrix(ModelingTransformMatrix, 4, 4, "ModelingTransformMatrix after reset");
}

// Multiplication Cross for 3-dimensionality
void Cross_Multi(float destination[], float a[], float b[]) {
	destination[0] = a[1] * b[2] - a[2] * b[1];
	destination[1] = a[2] * b[0] - a[0] * b[2];
	destination[2] = a[0] * b[1] - a[1] * b[0];
}
// Unitize a 3-dimensionality Vector
void UnitizeVector(float vector[]) {
	float diver = sqrt(vector[0] * vector[0] +
		vector[1] * vector[1] +
		vector[2] * vector[2]);
	for (int i = 0; i < 3; i++) {
		vector[i] /= diver;
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
// Two 4 * 4 Matrixs Multiplication
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

void scale(float sx, float sy, float sz) {
	float scaling_matrix[MTM_SIZE][MTM_SIZE]
		= { { sx, 0, 0, 0 },
			{ 0, sy, 0, 0 },
			{ 0, 0, sz, 0 },
			{ 0, 0, 0, 1 } };

	float new_matrix[MTM_SIZE][MTM_SIZE];

	Matrix_Multi_Matrix(scaling_matrix, ModelingTransformMatrix, new_matrix);

	update_Matrix(new_matrix, ModelingTransformMatrix);

	printMatrix(ModelingTransformMatrix, 4, 4, "ModelingTransformMatrix after scale");
}
void translate(float tx, float ty, float tz) {
	float translate_matrix[MTM_SIZE][MTM_SIZE]
		= { { 1, 0, 0, tx },
			{ 0, 1, 0, ty },
			{ 0, 0, 1, tz },
			{ 0, 0, 0, 1 } };

	float new_matrix[MTM_SIZE][MTM_SIZE];

	Matrix_Multi_Matrix(translate_matrix, ModelingTransformMatrix, new_matrix);

	update_Matrix(new_matrix, ModelingTransformMatrix);

	printMatrix(ModelingTransformMatrix, 4, 4, "ModelingTransformMatrix after transalte");
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
		printMatrix(ModelingTransformMatrix, 4, 4, "ModelingTransformMatrix after rotate");
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
		printMatrix(ModelingTransformMatrix, 4, 4, "ModelingTransformMatrix");
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
		printMatrix(ModelingTransformMatrix, 4, 4, "ModelingTransformMatrix");
	}
}

// Set up Eye_Transform_Matrix
void Eye_Transform(float PX, float PY, float PZ, float CX, float CY, float CZ, float Tilt) {
	float Eye_Translation_Matrix[MTM_SIZE][MTM_SIZE]
		= { { 1, 0, 0, -PX },
			{ 0, 1, 0, -PY },
			{ 0, 0, 1, -PZ },
			{ 0, 0, 0, 1 } };
	printMatrix(Eye_Translation_Matrix, 4, 4, "Eye_Translation_Matrix");

	float Eye_Mirror_Matrix[MTM_SIZE][MTM_SIZE]
		= { { -1, 0, 0, 0 },
			{ 0, 1, 0, 0 },
			{ 0, 0, 1, 0 },
			{ 0, 0, 0, 1 } };
	printMatrix(Eye_Mirror_Matrix, 4, 4, "Eye_Mirror_Matrix");

	float TiltDegree = Tilt * M_PI / 180.0;
	float Eye_Tilt_Matrix[MTM_SIZE][MTM_SIZE]
		= { { cos(TiltDegree), sin(TiltDegree), 0, 0 },
			{ -sin(TiltDegree), cos(TiltDegree), 0, 0 },
			{ 0, 0, 1, 0 },
			{ 0, 0, 0, 1 } };
	printMatrix(Eye_Tilt_Matrix, 4, 4, "Eye_Tilt_Matrix");

	float vectorZ[3] = { CX - PX, CY - PY, CZ - PZ };
	float top_vector[3] = { 0, 1, 0 };

	view_vactor[0] = vectorZ[0];
	view_vactor[1] = vectorZ[1];
	view_vactor[2] = vectorZ[2];
	float vector3[3];
	memcpy(vector3, vectorZ, sizeof(vector3));
	UnitizeVector(vector3);
	float vector1[3];
	Cross_Multi(vector1, top_vector, vectorZ);
	UnitizeVector(vector1);
	float vector2[3];
	Cross_Multi(vector2, vector3, vector1);
	UnitizeVector(vector2);

	float GRM[MTM_SIZE][MTM_SIZE]
		= { { vector1[0], vector1[1], vector1[2], 0 },
			{ vector2[0], vector2[1], vector2[2], 0 },
			{ vector3[0], vector3[1], vector3[2], 0 },
			{ 0, 0, 0, 1 } };
	printMatrix(GRM, 4, 4, "GRM");

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
	printMatrix(Eye_Matirx, 4, 4, "Eye_Matrix");
}
// Set up Project_Transform_Matrix
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

void drawFrame(int Frame[]) {
	drawLine(Frame[0], Frame[2], Frame[1], Frame[2], 1.0, 1.0, 1.0);
	num_line++;

	drawLine(Frame[1], Frame[2], Frame[1], Frame[3], 1.0, 1.0, 1.0);
	num_line++;

	drawLine(Frame[1], Frame[3], Frame[0], Frame[3], 1.0, 1.0, 1.0);
	num_line++;

	drawLine(Frame[0], Frame[3], Frame[0], Frame[2], 1.0, 1.0, 1.0);
	num_line++;
}
void object(string objectname) {
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

	int n;
	// read face one by one
	for (int i = 0; i < ASCModel[num_ASCModel].num_face; i++) {
		fin >> n;
		for (int j = 0; j < n; j++) {
			fin >> ASCModel[num_ASCModel].face[i][j];
		}
		if (n < MAX_MODEL) {
			for (int j = n; j < MAX_MODEL; j++)
				ASCModel[num_ASCModel].face[i][j] = -1;
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

	num_ASCModel++;
}
// Set the observer
void observer(float PX, float PY, float PZ, float CX, float CY, float CZ,
	float Tilt, float zNear, float zFar, float hFOV) {
	Eye_Transform(PX, PY, PZ, CX, CY, CZ, Tilt);

	Project_Transform(zNear, zFar, hFOV);

	WzNear = zNear;
	WzFar = zFar;
	WhFOV = hFOV;
}
// Set the view border
void viewport(float vxl, float vxr, float vyb, float vyt) {
	float Vlength = vxr - vxl;
	float Vheight = vyt - vyb;
	Project_Matrix[1][1] = Vlength / Vheight;// Set Ar
	printMatrix(Project_Matrix, 4, 4, "Project_Matrix");

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
	printMatrix(WVM, 4, 4, "WVM Matrix");
}

responseType clipProcess(deque<Point> source, const int code, const int faceIndex, const int faceSize) {
	deque<Point> target;
	float c1;
	float c2;
	bool isClip = false;
	for (int j = 0; j < source.size(); j++) {
		const Point point = source.at(j);
		const Point point1 = source.at((j + 1) % source.size());
		switch (code) {
		case WMinusX: {
			// w -x
			c1 = point.w - point.x;
			c2 = point1.w - point1.x;
			break;
		}
		case WPlusX: {
			// w +x
			c1 = point.w + point.x;
			c2 = point1.w + point1.x;
			break;
		}
		case WMinusY: {
			// w - y
			c1 = point.w - point.y;
			c2 = point1.w - point1.y;
			break;
		}
		case WPlusY: {
			// w + y
			c1 = point.w + point.y;
			c2 = point1.w + point1.y;
			break;
		}
		case WMinusZ: {
			// w -z
			c1 = point.w - point.z;
			c2 = point1.w - point1.z;
			break;
		}
		case Z0: {
			// z >= 0
			c1 = point.z;
			c2 = point1.z;
			break;
		}
		default: {
		}
		}

		try {
			if (c1 >= 0 && c2 >= 0) {
				Point temp;
				temp.x = point1.x;
				temp.y = point1.y;
				temp.z = point1.z;
				temp.w = point1.w;
				temp.vertexNumber = point1.vertexNumber;
				target.push_back(temp);
			}
			else if (c1 < 0 && c2 >= 0) {
				Point tempI;
				tempI.x = point.x + (c1 / (c1 - c2)) * (point1.x - point.x);
				tempI.y = point.y + (c1 / (c1 - c2)) * (point1.y - point.y);
				tempI.z = point.z + (c1 / (c1 - c2)) * (point1.z - point.z);
				tempI.w = point.w + (c1 / (c1 - c2)) * (point1.w - point.w);
				tempI.vertexNumber = Cur_ASCModel.num_vertex + 1;
				target.push_back(tempI);

				Point tempP;
				tempP.x = point1.x;
				tempP.y = point1.y;
				tempP.z = point1.z;
				tempP.w = point1.w;
				tempP.vertexNumber = point1.vertexNumber;
				target.push_back(tempP);
				isClip = true;
			}
			else if (c1 >= 0 && c2 < 0) {
				Point tempI;
				tempI.x = point.x + (c1 / (c1 - c2)) * (point1.x - point.x);
				tempI.y = point.y + (c1 / (c1 - c2)) * (point1.y - point.y);
				tempI.z = point.z + (c1 / (c1 - c2)) * (point1.z - point.z);
				tempI.w = point.w + (c1 / (c1 - c2)) * (point1.w - point.w);
				tempI.vertexNumber = Cur_ASCModel.num_vertex + 1;//point1.vertexNumber;
				target.push_back(tempI);
				isClip = true;
			}
		}
		catch (exception e) {
			cout << "check for this break point" << endl;
		}
	}
	responseType response;
	response.polygon = target;
	response.isClip = isClip;
	return response;
}
responseType clip(const int faceSize, const int faceIndex) {
	// inital data
	deque<Point> source;
	responseType response;
	bool isClip = false;
	for (int j = 0; j < faceSize; j++) {
		Point temp;
		temp.x = Cur_ASCModel.vertex[Cur_ASCModel.face[faceIndex][j] - 1][0];
		temp.y = Cur_ASCModel.vertex[Cur_ASCModel.face[faceIndex][j] - 1][1];
		temp.z = Cur_ASCModel.vertex[Cur_ASCModel.face[faceIndex][j] - 1][2];
		temp.w = Cur_ASCModel.vertex[Cur_ASCModel.face[faceIndex][j] - 1][3];
		temp.vertexNumber = Cur_ASCModel.face[faceIndex][j];
		source.push_back(temp);
	}

	// do clip
	//w -x
	response = clipProcess(source, WMinusX, faceIndex, faceSize);
	source.clear();
	source = response.polygon;
	isClip = isClip | response.isClip;

	//w + x
	response = clipProcess(source, WPlusX, faceIndex, faceSize);
	source.clear();
	source = response.polygon;
	isClip = isClip | response.isClip;

	//w - y
	response = clipProcess(source, WMinusY, faceIndex, faceSize);
	source.clear();
	source = response.polygon;
	isClip = isClip | response.isClip;

	//w + y
	response = clipProcess(source, WPlusY, faceIndex, faceSize);
	source.clear();
	source = response.polygon;
	isClip = isClip | response.isClip;

	//w - z
	response = clipProcess(source, WMinusZ, faceIndex, faceSize);
	source.clear();
	source = response.polygon;
	isClip = isClip | response.isClip;

	//Z >=0
	response = clipProcess(source, Z0, faceIndex, faceSize);
	isClip = isClip | response.isClip;
	response.isClip = isClip;

	return response;
}

void clearScreen() {
	glClear(GL_COLOR_BUFFER_BIT);
	glFlush();
}
void clearData() {
	num_line = 0;
	return;
}

// hiden surface  removal & trasfer to Screen Space and Draw
void postProcessRendering(const responseType response) {
	bool isClip = response.isClip;
	float polygon[MAX_MODEL][MTM_SIZE] = { 0 };
	int count = 0;
	for (int i = 0; i < response.polygon.size(); i++) {
		polygon[i][0] = response.polygon.at(i).x;

		polygon[i][1] = response.polygon.at(i).y;

		polygon[i][2] = response.polygon.at(i).z;

		polygon[i][3] = response.polygon.at(i).w;

		count++;
	}

	// perspective devide
	for (int i = 0; i < response.polygon.size(); i++) {
		for (int j = 0; j < MTM_SIZE; j++) {
			if (polygon[i][3] != 0) {
				polygon[i][j] /= polygon[i][3];
			}
		}
		polygon[i][2] = 1;//Fix to three-dimensional homogeneous
	}

	if (NoBackFace) {
		float x0, x1, y0, y1;
		x0 = polygon[1][0] - polygon[0][0];
		y0 = polygon[1][1] - polygon[0][1];
		x1 = polygon[2][0] - polygon[1][0];
		y1 = polygon[2][1] - polygon[1][1];
		if (-((x0 * y1) - (x1 * y0)) < 0) {
			return;
		}
	}

	// Convert Project-Space to Screen-Space
	for (int i = 0; i < response.polygon.size(); i++) {
		float Matrix[MTM_SIZE];
		for (int j = 0; j < MTM_SIZE; j++) {
			float multi_sum = 0;
			for (int k = 0; k < MTM_SIZE; k++) {
				multi_sum += polygon[i][k]
					* WVM[j][k];
			}
			Matrix[j] = multi_sum;
		}
		// Update(cover) the Object-Space
		for (int j = 0; j < MTM_SIZE; j++) {
			polygon[i][j] = Matrix[j];
		}
	}

	// Draw
	for (int i = 0; i < response.polygon.size(); i++) {
		int i_next = (i + 1) % response.polygon.size();
		int x0 = polygon[i][0];
		int x1 = polygon[i_next][0];
		int y0 = polygon[i][1];
		int y1 = polygon[i_next][1];
		drawLine(x0, y0, x1, y1, R, G, B);
		num_line++;
	}
	glFlush();
}

void display() {
	clearScreen();
	clearData();

	int Frame[] = { frame.vxl, frame.vxr,frame.vyb, frame.vyt };
	drawFrame(Frame);
	glFlush();

	float multi_sum = 0;
	float Matrix[MTM_SIZE];
	for (int index = 0; index < num_ASCModel; index++) {
		Cur_ASCModel = ASCModel[index];

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
				Cur_ASCModel.vertex[i][j] = Matrix[j];
			}
		}
		// do clip
		for (int i = 0; i < ASCModel[index].num_face; i++) {
			int pointNum = 0;
			for (int j = 0; j < MAX_MODEL; j++) {
				if (Cur_ASCModel.face[i][j] != -1) {
					pointNum++;
				}
				else {
					break;
				}
			}
			/// test postProcessRedering();
			/*deque<Point> source;
			for (int j = 0; j < pointNum; j++) {
				Point temp;
				temp.x = Cur_ASCModel.vertex[Cur_ASCModel.face[i][j] - 1][0];
				temp.y = Cur_ASCModel.vertex[Cur_ASCModel.face[i][j] - 1][1];
				temp.z = Cur_ASCModel.vertex[Cur_ASCModel.face[i][j] - 1][2];
				temp.w = Cur_ASCModel.vertex[Cur_ASCModel.face[i][j] - 1][3];
				temp.vertexNumber = Cur_ASCModel.face[i][j];
				source.push_back(temp);
			}
			responseType response;
			response.polygon = source;
			response.isClip = false;
			postProcessRendering(response);*/
			postProcessRendering(clip(pointNum, i));
			//Sleep(15);
		}
	}
	cout << "\tObjects display successfully" << endl;
}

void readFile() {
	float sx, sy, sz;
	float Xdegree, Ydegree, Zdegree;
	float tx, ty, tz;
	float vxl, vxr, vyt, vyb;
	float PX, PY, PZ, CX, CY, CZ, Tilt, zNear, zFar, hFOV;
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
			object(objectname);
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
			system("pause");
		}

		else if (command == "end") {
			fin.close();
			system("pause");
			exit(0);
			return;
		}
		else if (command == "#") {
			getline(fin, comment);
			cout << "\t" << comment << endl;
		}
		else if (command == "nobackfaces") {
			NoBackFace = true;
			cout << "noBackFaces mode on" << endl;
		}
		cout << endl;
	}
}

void displayFunc() {
	readFile();
}
int main(int argc, char **argv) {
	// init GLUT and create Window
	string filename = argc >= 2 ? argv[1] : "Lab3B.in";
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
	glutCreateWindow("Assignment3");
	init();
	glutDisplayFunc(displayFunc);
	glutMainLoop();
	return 0;
}