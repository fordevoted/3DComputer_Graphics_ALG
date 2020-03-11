#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <GL/glut.h>
#include <deque>

using namespace std;

#define MATRIX_SIZE 3
#define SQUARE_SIZE 4

struct Point {
	float x, y;
	float z;
};

struct square_struct {
	Point right_up;
	Point left_up;
	Point left_down;
	Point right_down;
}SquareStore[30];
typedef struct square_struct squareStruct;
struct triangle_struct {
	Point first;
	Point second;
	Point third;
}TriangleStore[30];
typedef struct triangle_struct triangleStruct;
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

int width = 800;
int height = 600;
int num_square = 0;
int num_triangle = 0;
int num_line = 0;

int WVM[MATRIX_SIZE][MATRIX_SIZE] = { {15, 0, 100 }, { 0, 15, 100 }, { 0, 0, 1 } };
int Frame[4];
string filename;

// transform matrix
float current_matrix[MATRIX_SIZE][MATRIX_SIZE];

void printAll() {
	cout << "num_square: " << num_square << " num_triangle: " << num_triangle << " num_line: " << num_line << endl;
}
int sign(int a)
{
	return a > 0 ? 1 : (a < 0 ? -1 : 0);
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
	cout << message + " is" << endl;
	for (int i = 0; i < row; i++) {
		cout << "\t[ ";
		for (int j = 0; j < col; j++) {
			cout << matrix[i][j] << ", ";
		}
		cout << " ]" << endl;
	}
}
// 3 * 3 matrix multiplication
void matrixMultiplication(float a[][MATRIX_SIZE], float b[][MATRIX_SIZE], float c[][MATRIX_SIZE]) {
	float multi_sum = 0;
	for (int i = 0; i < MATRIX_SIZE; i++) {
		for (int j = 0; j < MATRIX_SIZE; j++) {
			multi_sum = 0;
			for (int k = 0; k < MATRIX_SIZE; k++) {
				multi_sum += b[k][j] * a[i][k];
			}
			c[i][j] = multi_sum;
		}
	}
}
void update_CurMatrix(float resource[][MATRIX_SIZE], float destination[][MATRIX_SIZE]) {
	for (int i = 0; i < MATRIX_SIZE; i++) {
		for (int j = 0; j < MATRIX_SIZE; j++) {
			destination[i][j] = resource[i][j];
		}
	}
}

void scale(float sx, float sy) {
	float scaling_matrix[MATRIX_SIZE][MATRIX_SIZE]
		= { { sx, 0, 0 },
			{ 0, sy, 0 },
			{ 0, 0, 1  } };

	float new_matrix[MATRIX_SIZE][MATRIX_SIZE];

	matrixMultiplication(scaling_matrix, current_matrix, new_matrix);

	update_CurMatrix(new_matrix, current_matrix);

	printMatrix(current_matrix, MATRIX_SIZE, MATRIX_SIZE, "current_matrix");
}
void rotate(float degree) {
	float rad = degree * M_PI / 180.0;
	float rotation_matrix[MATRIX_SIZE][MATRIX_SIZE]
		= { { cos(rad), -sin(rad), 0 },
			{ sin(rad), cos(rad),  0 },
			{ 0,        0,         1 } };

	//printMatrix(rotation_matrix, 3, 3, "rotation_matrix");

	float new_matrix[MATRIX_SIZE][MATRIX_SIZE];

	matrixMultiplication(rotation_matrix, current_matrix, new_matrix);

	update_CurMatrix(new_matrix, current_matrix);

	printMatrix(current_matrix, MATRIX_SIZE, MATRIX_SIZE, "current_matrix");
}
void translate(float tx, float ty) {
	float translate_matrix[MATRIX_SIZE][MATRIX_SIZE]
		= { { 1, 0, tx },
			{ 0, 1, ty },
			{ 0, 0, 1  } };

	float new_matrix[MATRIX_SIZE][MATRIX_SIZE];

	matrixMultiplication(translate_matrix, current_matrix, new_matrix);

	update_CurMatrix(new_matrix, current_matrix);

	printMatrix(current_matrix, MATRIX_SIZE, MATRIX_SIZE, "current_matrix");
}
void reset() {
	for (int i = 0; i < MATRIX_SIZE; i++) {
		for (int j = 0; j < MATRIX_SIZE; j++) {
			if (i == j)
				current_matrix[i][j] = 1.0;
			else
				current_matrix[i][j] = 0;
		}
	}
	printMatrix(current_matrix, 3, 3, "current_matrix");
}

void AddtoSquareStore(float matrix[][MATRIX_SIZE]) {
	SquareStore[num_square].right_up.x = matrix[0][0];
	SquareStore[num_square].right_up.y = matrix[0][1];
	SquareStore[num_square].right_up.z = matrix[0][2];

	SquareStore[num_square].left_up.x = matrix[1][0];
	SquareStore[num_square].left_up.y = matrix[1][1];
	SquareStore[num_square].left_up.z = matrix[1][2];

	SquareStore[num_square].left_down.x = matrix[2][0];
	SquareStore[num_square].left_down.y = matrix[2][1];
	SquareStore[num_square].left_down.z = matrix[2][2];

	SquareStore[num_square].right_down.x = matrix[3][0];
	SquareStore[num_square].right_down.y = matrix[3][1];
	SquareStore[num_square].right_down.z = matrix[3][2];
}
void square() {
	float orginal_square[4][MATRIX_SIZE]
		= { { 1,  1,  1 },
			{ -1, 1,  1 },
			{ -1, -1, 1 },
			{ 1,  -1, 1 } };

	float WorldCoord_square[4][MATRIX_SIZE];

	float multi_sum = 0;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < MATRIX_SIZE; j++) {
			multi_sum = 0;
			for (int k = 0; k < MATRIX_SIZE; k++) {
				multi_sum += orginal_square[i][k] * current_matrix[j][k];
			}
			WorldCoord_square[i][j] = multi_sum;
		}
	}

	AddtoSquareStore(WorldCoord_square);
	num_square++;

	printMatrix(WorldCoord_square, 4, MATRIX_SIZE, "WorldCoord_square");
}
void AddtoTriangleStore(float matrix[][MATRIX_SIZE]) {
	TriangleStore[num_triangle].first.x = matrix[0][0];
	TriangleStore[num_triangle].first.y = matrix[0][1];
	TriangleStore[num_triangle].first.z = matrix[0][2];

	TriangleStore[num_triangle].second.x = matrix[1][0];
	TriangleStore[num_triangle].second.y = matrix[1][1];
	TriangleStore[num_triangle].second.z = matrix[1][2];

	TriangleStore[num_triangle].third.x = matrix[2][0];
	TriangleStore[num_triangle].third.y = matrix[2][1];
	TriangleStore[num_triangle].third.z = matrix[2][2];
}
void triangle() {
	float org_triangle[3][MATRIX_SIZE]
		= { { 0,  1,  1 },
			{ -1, -1, 1 },
			{ 1,  -1, 1 } };

	float WorldCoord_square[3][MATRIX_SIZE];

	float multi_sum = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < MATRIX_SIZE; j++) {
			multi_sum = 0;
			for (int k = 0; k < MATRIX_SIZE; k++) {
				multi_sum += org_triangle[i][k] * current_matrix[j][k];
			}
			WorldCoord_square[i][j] = multi_sum;
		}
	}

	AddtoTriangleStore(WorldCoord_square);
	num_triangle++;

	printMatrix(WorldCoord_square, 3, MATRIX_SIZE, "WorldCoord_triangle");
}

void clearData() {
	num_line = 0;
	num_square = 0;
	num_triangle = 0;
	return;
}
void clearScreen() {
	glClear(GL_COLOR_BUFFER_BIT);
	glFlush();
}

void update_WVM(int sx, int sy, int tx, int ty) {
	WVM[0][0] = sx;
	WVM[1][1] = sy;
	WVM[0][2] = tx;
	WVM[1][2] = ty;
	WVM[2][2] = 1;
	WVM[0][1] = WVM[1][0] = WVM[2][0] = WVM[2][1] = 0;
}
void GetFromSquareStore(int index, float destination[][MATRIX_SIZE]) {
	destination[0][0] = SquareStore[index].right_up.x;
	destination[0][1] = SquareStore[index].right_up.y;
	destination[0][2] = SquareStore[index].right_up.z;

	destination[1][0] = SquareStore[index].left_up.x;
	destination[1][1] = SquareStore[index].left_up.y;
	destination[1][2] = SquareStore[index].left_up.z;

	destination[2][0] = SquareStore[index].left_down.x;
	destination[2][1] = SquareStore[index].left_down.y;
	destination[2][2] = SquareStore[index].left_down.z;

	destination[3][0] = SquareStore[index].right_down.x;
	destination[3][1] = SquareStore[index].right_down.y;
	destination[3][2] = SquareStore[index].right_down.z;
}
void GetFromTriangleStore(int index, float destination[][MATRIX_SIZE]) {
	destination[0][0] = TriangleStore[index].first.x;
	destination[0][1] = TriangleStore[index].first.y;
	destination[0][2] = TriangleStore[index].first.z;

	destination[1][0] = TriangleStore[index].second.x;
	destination[1][1] = TriangleStore[index].second.y;
	destination[1][2] = TriangleStore[index].second.z;

	destination[2][0] = TriangleStore[index].third.x;
	destination[2][1] = TriangleStore[index].third.y;
	destination[2][2] = TriangleStore[index].third.z;
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
void DrawFrame() {
	drawLine(Frame[0], Frame[2], Frame[1], Frame[2], 1.0, 1.0, 1.0);
	num_line++;

	drawLine(Frame[1], Frame[2], Frame[1], Frame[3], 1.0, 1.0, 1.0);
	num_line++;

	drawLine(Frame[1], Frame[3], Frame[0], Frame[3], 1.0, 1.0, 1.0);
	num_line++;

	drawLine(Frame[0], Frame[3], Frame[0], Frame[2], 1.0, 1.0, 1.0);
	num_line++;
}

responseType clip(int matrix[][MATRIX_SIZE], int PointNum, struct frame_struct Frame) {
	deque<Point> source;
	deque<Point> target;
	bool isClip = false;

	// inital data
	for (int i = 0; i < PointNum; i++) {
		Point temp;
		temp.x = matrix[i][0];
		temp.y = matrix[i][1];
		source.push_back(temp);
	}

	// do clip
	float b;
	float a;

	//left
	for (int j = 0; j < source.size(); j++) {
		Point point = source.at(j);
		Point point1 = source.at((j + 1) % source.size());
		if (point.x != point1.x) {
			b = (point.x * point1.y - point.y * point1.x)
				/ (point.x - point1.x);
			a = (point.y - point1.y) / (point.x - point1.x);
		}
		else {
			b = point.x;
			a = INFINITY;
		}
		if (point.x >= Frame.vxl
			&&point1.x >= Frame.vxl) {
			Point temp;
			temp.x = point1.x;
			temp.y = point1.y;
			target.push_back(temp);
		}
		else if (point.x < Frame.vxl
			&&point1.x > Frame.vxl) {
			Point tempI;
			tempI.x = Frame.vxl;
			tempI.y = a * Frame.vxl + b;
			target.push_back(tempI);
			Point tempP;
			tempP.x = point1.x;
			tempP.y = point1.y;
			target.push_back(tempP);
			isClip = true;
		}
		else if (point.x > Frame.vxl
			&&point1.x < Frame.vxl) {
			Point temp;
			temp.x = Frame.vxl;
			temp.y = a * Frame.vxl + b;
			target.push_back(temp);
			isClip = true;
		}
	}

	source.clear();
	for (int j = 0; j < target.size(); j++) {
		source.push_back(target.at(j));
	}
	target.clear();

	//bot
	for (int j = 0; j < source.size(); j++) {
		Point point = source.at(j);
		Point point1 = source.at((j + 1) % source.size());
		;		if (point.x != point1.x) {
			b = (point.x * point1.y - point.y * point1.x) /
				(point.x - point1.x);
			a = (point.y - point1.y) / (point.x - point1.x);
		}
		else {
			b = point.x;
			a = INFINITY;
		}
		if (point.y >= Frame.vyb
			&&point1.y >= Frame.vyb) {
			Point temp;
			temp.x = point1.x;
			temp.y = point1.y;
			target.push_back(temp);
		}
		else if (point.y < Frame.vyb
			&&point1.y > Frame.vyb) {
			Point tempI;
			if (a == INFINITY || a == 0) {
				tempI.x = point1.x;
			}
			else {
				tempI.x = (Frame.vyb - b) / a;
			}
			tempI.y = Frame.vyb;
			target.push_back(tempI);
			Point tempP;
			tempP.x = point1.x;
			tempP.y = point1.y;
			target.push_back(tempP);
			isClip = true;
		}
		else if (point.y > Frame.vyb
			&&point1.y < Frame.vyb) {
			Point temp;
			if (a == INFINITY || a == 0) {
				temp.x = point1.x;
			}
			else {
				temp.x = (Frame.vyb - b) / a;
			}
			temp.y = Frame.vyb;
			target.push_back(temp);
			isClip = true;
		}
	}

	source.clear();
	for (int j = 0; j < target.size(); j++) {
		source.push_back(target.at(j));
	}
	target.clear();

	//right
	for (int j = 0; j < source.size(); j++) {
		Point point = source.at(j);
		Point point1 = source.at((j + 1) % source.size());
		;		if (point.x != point1.x) {
			b = (point.x * point1.y - point.y * point1.x) /
				(point.x - point1.x);
			a = (point.y - point1.y) / (point.x - point1.x);
		}
		else {
			b = point.x;
			a = INFINITY;
		}
		if (point.x <= Frame.vxr
			&&point1.x <= Frame.vxr) {
			Point temp;
			temp.x = point1.x;
			temp.y = point1.y;
			target.push_back(temp);
		}
		else if (point.x > Frame.vxr
			&&point1.x < Frame.vxr) {
			Point tempI;
			tempI.x = Frame.vxr;
			tempI.y = a * Frame.vxr + b;
			target.push_back(tempI);
			Point tempP;
			tempP.x = point1.x;
			tempP.y = point1.y;
			target.push_back(tempP);
			isClip = true;
		}
		else if (point.x < Frame.vxr
			&&point1.x > Frame.vxr) {
			Point temp;
			temp.x = Frame.vxr;
			temp.y = a * Frame.vxr + b;
			target.push_back(temp);
			isClip = true;
		}
	}

	source.clear();
	for (int j = 0; j < target.size(); j++) {
		source.push_back(target.at(j));
	}
	target.clear();

	//top
	for (int j = 0; j < source.size(); j++) {
		Point point = source.at(j);
		Point point1 = source.at((j + 1) % source.size());
		if (point.x != point1.x) {
			b = (point.x * point1.y - point.y * point1.x) /
				(point.x - point1.x);
			a = (point.y - point1.y) / (point.x - point1.x);
		}
		else {
			b = point.x;
			a = INFINITY;
		}
		if (point.y <= Frame.vyt
			&&point1.y <= Frame.vyt) {
			Point temp;
			temp.x = point1.x;
			temp.y = point1.y;
			target.push_back(temp);
		}
		else if (point.y > Frame.vyt
			&&point1.y < Frame.vyt) {
			Point tempI;
			if (a == INFINITY || a == 0) {
				tempI.x = point1.x;
			}
			else {
				tempI.x = (Frame.vyt - b) / a;
			}
			tempI.y = Frame.vyt;
			target.push_back(tempI);
			Point tempP;
			tempP.x = point1.x;
			tempP.y = point1.y;
			target.push_back(tempP);
			isClip = true;
		}
		else if (point.y < Frame.vyt
			&&point1.y > Frame.vyt) {
			Point temp;
			if (a == INFINITY || a == 0) {
				temp.x = point1.x;
			}
			else {
				temp.x = (Frame.vyt - b) / a;
			}
			temp.y = Frame.vyt;
			target.push_back(temp);
			isClip = true;
		}
	}

	responseType response;
	response.isClip = isClip;
	response.polygon = target;
	return response;
}
void DrawSquareDirect(int matrix[][MATRIX_SIZE], struct frame_struct Frame) {
	responseType response = clip(matrix, 4, Frame);
	for (int i = 0; i < response.polygon.size(); i++) {
		if (i == response.polygon.size() - 1) {
			drawLine(response.polygon.at(i).x, response.polygon.at(i).y,
				response.polygon.at(0).x, response.polygon.at(0).y,
				1.0, 0.5, 0.8);
			break;
		}
		drawLine(response.polygon.at(i).x, response.polygon.at(i).y,
			response.polygon.at(i + 1).x, response.polygon.at(i + 1).y,
			1.0, 0.5, 0.8);
		num_line++;
	}
}
void DrawTriangleDirect(int matrix[][MATRIX_SIZE], struct frame_struct Frame) {
	responseType response = clip(matrix, 3, Frame);
	for (int i = 0; i < response.polygon.size(); i++) {
		if (i == response.polygon.size() - 1) {
			drawLine(response.polygon.at(i).x, response.polygon.at(i).y,
				response.polygon.at(0).x, response.polygon.at(0).y,
				1.0, 0.5, 0.8);
			num_line++;
			break;
		}
		drawLine(response.polygon.at(i).x, response.polygon.at(i).y,
			response.polygon.at(i + 1).x, response.polygon.at(i + 1).y,
			1.0, 0.5, 0.8);
		num_line++;
	}
}

void view(float wxl, float wxr, float wyb, float wyt, float vxl, float vxr, float vyb, float vyt) {
	int scaling_x = (int)((vxr - vxl) / (wxr - wxl));
	int scaling_y = (int)((vyt - vyb) / (wyt - wyb));
	int shift_x = (int)(vxl - scaling_x * wxl);
	int shift_y = (int)(vyb - scaling_y * wyb);
	update_WVM(scaling_x, scaling_y, shift_x, shift_y);

	printMatrix(WVM, 3, 3, "WVM");

	float WorldCoord_matrix[4][MATRIX_SIZE];
	int ScreenCoord_matrix[4][MATRIX_SIZE];

	//Set Frame
	struct frame_struct frame;
	Frame[0] = frame.vxl = (int)vxl;
	Frame[1] = frame.vxr = (int)vxr;
	Frame[2] = frame.vyb = (int)vyb;
	Frame[3] = frame.vyt = (int)vyt;

	// Draw Frame
	DrawFrame();

	//Draw Square
	float multi_sum = 0;
	int count = 0;
	while (count < num_square) {
		GetFromSquareStore(count, WorldCoord_matrix);
		printMatrix(WorldCoord_matrix, 4, 3, "WorldCoord_matrix");
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < MATRIX_SIZE; j++) {
				multi_sum = 0;
				for (int k = 0; k < MATRIX_SIZE; k++) {
					multi_sum += WorldCoord_matrix[i][k] * WVM[j][k];
				}
				ScreenCoord_matrix[i][j] = (int)multi_sum;
			}
		}
		printMatrix(ScreenCoord_matrix, 4, 3, "ScreenCoord_matrix");
		DrawSquareDirect(ScreenCoord_matrix, frame);
		//void DrawSquareByindex(count);
		count++;
	}

	//Draw Triangle
	count = 0;
	while (count < num_triangle) {
		GetFromTriangleStore(count, WorldCoord_matrix);
		printMatrix(WorldCoord_matrix, 3, 3, "WorldCoord_matrix");
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < MATRIX_SIZE; j++) {
				multi_sum = 0;
				for (int k = 0; k < MATRIX_SIZE; k++) {
					multi_sum += WorldCoord_matrix[i][k] * WVM[j][k];
				}
				ScreenCoord_matrix[i][j] = (int)multi_sum;
			}
		}
		printMatrix(ScreenCoord_matrix, 3, 3, "ScreenCoord_matrix");
		DrawTriangleDirect(ScreenCoord_matrix, frame);
		count++;
	}

	glFlush();
	return;
}

void ReadFile(string filename) {
	/*#if READ_MULTI_FILE
		ifstream fin(file_list[cur_file]);
	#else
		ifstream fin("lab2A.in");
	#endif*/
	ifstream fin(filename);
	if (fin.is_open()) {
		cout << "open the file successfully" << endl;
	}
	float sx, sy, degree, tx, ty, wxl, wxr, wyb, wyt, vxl, vxr, vyt, vyb;
	string command, comment;
	while (!fin.eof()) {
		fin >> command;
		if (command == "scale") {
			fin >> sx;
			fin >> sy;
			scale(sx, sy);
			cout << command << endl;
		}
		else if (command == "rotate") {
			fin >> degree;
			rotate(degree);
			cout << command << endl;
		}
		else if (command == "translate") {
			fin >> tx;
			fin >> ty;
			translate(tx, ty);
			cout << command << endl;
		}
		else if (command == "reset") {
			reset();
			cout << command << endl;
		}
		else if (command == "square") {
			square();
			cout << command << endl;
		}
		else if (command == "triangle") {
			triangle();
			cout << command << endl;
		}
		else if (command == "view") {
			fin >> wxl >> wxr >> wyb >> wyt >> vxl >> vxr >> vyb >> vyt;
			view(wxl, wxr, wyb, wyt, vxl, vxr, vyb, vyt);
			printAll();
			cout << command << endl;
			system("pause");
		}
		else if (command == "clearData") {
			clearData();
			cout << command << endl;
		}
		else if (command == "clearScreen") {
			clearScreen();
			cout << command << endl;
		}
		else if (command == "end") {
			exit(0);
			return;
		}
		else if (command == "#") {
			cout << command << endl;
			getline(fin, comment);
		}
	}
}
void display() {
	cout << "in display" << endl;
	ReadFile(filename);
}

int main(int argc, char **argv) {
	// init GLUT and create Window
	if (argc >= 2) {
		filename = argv[1];
	}
	else {
		filename = "lab2E.in";
	}
	system("pause");
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(width, height);
	glutCreateWindow("Assignment1");
	init();
	glutDisplayFunc(display);
	/*int m[4][3] = { {70,36,1},{65,30,1},{70,24,1},{76,30,1} };
	int m2[4][3] = { {10,12.8,1},{7.17,10,1},{10,7.17,1},{12.8,10,1} };

	DrawSquareDirect(m);
	glFlush();*/
	glutMainLoop();
	return 0;
}