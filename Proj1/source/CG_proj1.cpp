#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

using Matrix3d = std::array<std::array<double, 3>, 3>;
using Vector3d = std::array<double, 3>;
using Coord = std::array<int, 2>;
constexpr int POINTSIZE{ 1 };	// the number of pixels a point maps to
constexpr int WIDTH{ 800 };
constexpr int HEIGHT{ 600 };

std::string file_path;
std::ifstream infile;
std::stringstream ss;

constexpr Matrix3d I{ {
						{1, 0, 0},
						{0, 1, 0},
						{0, 0, 1}
					 } };
constexpr double PI{ 3.14159265358979323846 };
constexpr double PIdiv180{ PI / 180 };

struct ops {					// ops stands for operations
	unsigned char op;
	std::vector<Coord> mouseCoords;
};
std::vector<Coord> mouseCoords;	// cleaned up after every draw and upon key press
std::vector<ops>   drawLog;		// records all drawing activities
std::vector<Coord> polyCoords;	// used for drawing polygons exclusively

namespace cexp {
	constexpr void swap(double& a, double& b) {
		const double temp{ a };
		a = b;
		b = temp;
	}
}

// transpose a matrix
// Matrix3d ⟼ Matrix3d
constexpr Matrix3d transpose(const Matrix3d& A) {
	Matrix3d M{ A };
	cexp::swap(M[0][1], M[1][0]);
	cexp::swap(M[0][2], M[2][0]);
	cexp::swap(M[1][2], M[2][1]);
	return M;
}

// standard inner product
// Vector3d X Vector3d ⟼ double
constexpr double operator*(const Vector3d& lhs, const Vector3d& rhs) {
	return (lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2]);
}

// Matrix3d X Matrix3d ⟼ Matrix3d
constexpr Matrix3d operator*(const Matrix3d& lhs, const Matrix3d& rhs) {
	Matrix3d product{ {{0}} };
	Matrix3d rhst{ transpose(rhs) };
	for (size_t i = 0; i < 3; ++i) {
		for (size_t j = 0; j < 3; ++j) {
			product[i][j] = lhs[i] * rhst[j];
		}
	}
	return product;
}

// Matrix3d X Vector3d ⟼ Vector3d
constexpr Vector3d operator*(const Matrix3d& lhs, const Vector3d& rhs) {
	return { lhs[0] * rhs,
			lhs[1] * rhs,
			lhs[2] * rhs };
}

// double X Matrix3d ⟼ Vector3d
constexpr Vector3d operator*(const double lhs, const Vector3d& rhs) {
	return { lhs * rhs[0],
			lhs * rhs[1],
			lhs * rhs[2] };
}

constexpr std::ostream& operator<<(std::ostream& out, const Vector3d& v) {
	for (const auto& a : v) {
		out << "[ " << a << " ]" << '\n';
	}
	return out;
}

constexpr std::ostream& operator<<(std::ostream& out, const Matrix3d& m) {
	for (const auto& row : m) {
		out << "[ ";
		for (const auto& scalar : row) {
			out << scalar << " ";
		}
		out << "]\n";
	}
	return out;
}

// return a scaling matrix
constexpr Matrix3d getSM(const double Sx, const double Sy) {
	return { {
		{Sx, 0, 0},
		{0, Sy, 0},
		{0,  0, 1}
	} };
}

// return a translation matrix
constexpr Matrix3d getTM(const double dx, const double dy) {
	return { {
		{1, 0, dx},
		{0, 1, dy},
		{0, 0,  1}
	} };
}

// return a rotation matrix
inline Matrix3d getRM(const double radian) {
	const double cosine{ std::cos(radian) };
	const double sine{ std::sin(radian) };
	return { {
		{cosine,   -sine, 0},
		{  sine,  cosine, 0},
		{     0,       0, 1}
	} };
}

// simple hash function
constexpr size_t operator""_hash(const char* s, const size_t count) {
	return *s ^ count;
}

void drawALine(Coord endpoint1, Coord endpoint2) {
	// unpack data
	auto [x1, y1] { endpoint1 };		// structured binding,requires C++17. didn't know u could do this, cool.
	auto [x2, y2] { endpoint2 };
	// drawing logic begins
	if (x2 < x1) {	                    // let (x1, y1) be the point on the left
		std::swap(x1, x2);
		std::swap(y1, y2);
	}
	glPointSize(POINTSIZE);
	glBegin(GL_POINTS);
	glVertex2i(x1, y1);			// draw the first point no matter wut

	//preprocessing
	bool negativeSlope{ y1 > y2 };
	if (negativeSlope) {		// see if the slope is negative
		y2 += 2 * (y1 - y2);	// mirror the line with respect to y = y1 for now, and voila, a positive slope line!
	}							// we draw this line as if the slope was positive in our head and draw it "upside down" on the screen in the while loop
	bool mGreaterThanOne{ (y2 - y1) > (x2 - x1) };
	if (mGreaterThanOne) {		// slope greater than 1, swap x and y,  mirror the line with respect to y = x for now, mirror it again when drawing
		std::swap(x1, y1);
		std::swap(x2, y2);
	}
	int x{ x1 };
	int y{ y1 };
	int a{ y2 - y1 };
	int b{ x1 - x2 };
	int d{ 2 * a + b };
	while (x < x2) {		// draw from left to right (recall that we make x2 be always on the right)
		if (d <= 0) {		// choose E
			d += 2 * a;
			if (mGreaterThanOne) {							// sort of "mirror" the negative slope line back to where it's supposed to be
				glVertex2i(y, !negativeSlope ? (++x) : (2 * x1 - ++x)); // slope > 1 y is actually x and vice versa
			}
			else {
				glVertex2i(++x, !negativeSlope ? (y) : (2 * y1 - y));
			}
		}
		else {				// choose NE
			d += 2 * (a + b);
			if (mGreaterThanOne) {
				glVertex2i(++y, !negativeSlope ? (++x) : (2 * x1 - ++x));
			}
			else {
				glVertex2i(++x, !negativeSlope ? (++y) : (2 * y1 - ++y));
			}
		}
	}
	glEnd();
}

inline void drawPolyAtOnce(const std::vector<Coord>& V) {
	drawALine(V.front(), V[1]);
	drawALine(V.front(), V.back());
	for (auto it = V.cbegin() + 1; it != V.cend(); it++) {
		drawALine(*(it - 1), *it);
	}
}

void handleLine() {
}

void displayFunc() {
	std::string line, str;
	Matrix3d T{ I };
	double Sx, Sy,
		degree, radian,
		dx, dy,
		wxl, wxr, wyb, wyt, vxl, vxr, vyb, vyt;
	while (std::getline(infile, line)) {
		ss << line;
		ss >> str;
		switch (operator""_hash(str.c_str(), str.size())) {
		case "scale"_hash:          // apply scaling
			ss >> Sx >> Sy;
			std::cout << "\nScale, Sx=" << Sx << ", Sy=" << Sy
				<< "\nbefore updating\n"
				<< T
				<< "after updating\n";
			T = getSM(Sx, Sy) * T;
			std::cout << T;
			break;
		case "rotate"_hash:         // apply rotation
			ss >> degree;
			radian = degree * PIdiv180;  // convert degree to radian
			std::cout << "\nRotate, degree=" << degree
				<< ", radian=" << radian
				<< "\nbefore updating\n"
				<< T
				<< "after updating\n";
			T = getRM(radian) * T;
			std::cout << T;
			break;
		case "translate"_hash:      // apply translation
			ss >> dx >> dy;
			std::cout << "\nTranslate, dx=" << dx << ", dy=" << dy
				<< "\nbefore updating\n"
				<< T
				<< "after updating\n";
			T = getTM(dx, dy) * T;
			std::cout << T;
			break;
		case "square"_hash:         // draw a square
			std::cout << "\nOne square created\n";
			drawALine({ 0, 0 }, { 600, 600 });
			break;
		case "triangle"_hash:       // draw a triangle
			std::cout << "\nOne triangle created\n";
			drawALine({ 600, 0 }, { 600, 600 });
			break;
		case "view"_hash:           // create a view (map to the screen)
			// To-do: print some info and draw on glut windows
			// bottom left corner is the origin

			glFlush();
			std::cout << "(press enter to continue...)";
			getchar();  // system("pause"); on Windows
			break;
		case "clearData"_hash:      // clear all recorded data, zero out all stats
			std::cout << "\nClear all data.\n";
			break;
		case "clearScreen"_hash:    // clear glut window
			std::cout << "\nClear screen.\n";
			break;
		case "end"_hash:            // terminate the process
			exit(EXIT_SUCCESS);
		case "reset"_hash:          // set transformation matrix to identity
			std::cout << "\nReset transformation.\n";
			T = I;
			break;
		case ""_hash:
		case "#"_hash:
			break;
		default:
			exit(EXIT_FAILURE);
		}
		while (ss >> str);
		str.clear();
		ss.clear();
	}
}

int main(int argc, char* argv[]) {
	std::cout << "(press enter to continue...)";
	getchar();

	file_path = ((argc == 2) ? argv[1] : "lab2A.in");
	infile.open(file_path);

	// init GLUT and create Window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(WIDTH, HEIGHT);
	glutCreateWindow("Your First GLUT Window!");
	// displayFunc is called whenever there is a need to redisplay the
	// window, e.g. when the window is exposed from under another window or
	// when the window is de-iconified
	glutDisplayFunc(displayFunc);
	// set bargckground color
	glClearColor(0.0, 0.0, 0.0, 0.0); // set the bargckground
	glClear(GL_COLOR_BUFFER_BIT); // clear the buffer
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, WIDTH, 0.0, HEIGHT);
	glFlush();
	// enter GLUT event processing cycle
	glutMainLoop();
	std::cout << "end\n";
}