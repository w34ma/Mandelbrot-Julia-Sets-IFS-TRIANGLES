//
// IMPORTANT: fill in things below BEFORE submitting your solution
//
//-------------------------------------------------------------------------------------
// <course_id> CSIS0271_COMP3271 </course_id>
//
// <course_name> Computer Graphics <course_name>
//
// <student_name> Weicong Ma </student_name>
//
// <student_id> 3035211244 </student_id>
//
// <submission_date> 2015-02-14 </submission_date>
//
// <notes_made_to_us> fix points in  triangle_list</notes_made_to_us>
//
// <assignment_description> refer to PDF file </assignment_description>
//-------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------
//  cg_assignment_1.cpp
//
//  cg_assignment_1_template - template for COMP3271 Computer Graphics assignment 1
//
//  Instructor: Jack M. Wang
//
//  TA: Kan WU, ulmonkey1987@gmail.com, kwu@cs.hku.hk
//
//  Department of Computer Science, The University of Hong Kong
//-------------------------------------------------------------------------------------

#if defined(_WIN32)

    #include <windows.h>
    #include <GL/gl.h>
    #include <GL/glu.h>
    #include <GL/glut.h>

#elif defined(__linux)

    #include <GL/gl.h>
    #include <GL/glu.h>
    #include <GL/glut.h>

#elif defined(__APPLE__)

    #include <OpenGL/gl.h>
    #include <OpenGL/glu.h>
    #include <GLUT/glut.h>

#endif

#include <stdint.h>

#include <iostream>
using std::cout;
using std::endl;

#include <complex>
using std::complex;

#include <vector>
using std::vector;

#include <cmath>
//using std::min;
//using std::max;

//-------------------------------------------------------------------------------------
// data type definition

// assignment type
// change this will change the current displayed things:
enum enumAssignmentType
{
    MANDABROT_SET = 0, JULIA_SET, IFS_TRIANGLES
};

// 3-component double-precision color (r, g, b: [0.0, 1.0])
class color3f
{
public:
    
    color3f(void): r_(0.0), g_(0.0), b_(0.0) {};
    
    color3f(const double r, const double g, const double b)
            : r_(r), g_(g), b_(b) {};
    
    double r_, g_, b_;
};

// 2d double-precision points
class point2f
{
public:
    
    point2f(void): x_(0), y_(0) {};
    
    point2f(const double x, const double y): x_(x), y_(y) {};
    
    point2f(const point2f& p): x_(p.x_), y_(p.y_) {};
    
    double x_, y_;
};

// triangle type
class triangle_t
{
public:
    
    triangle_t(void) {};
    
    triangle_t(const point2f& p1, const point2f& p2, const point2f& p3)
                : p1_(p1), p2_(p2), p3_(p3) {};
    
    triangle_t(const double x1, const double y1,
               const double x2, const double y2,
               const double x3, const double y3)
                : p1_(x1, y1), p2_(x2, y2), p3_(x3, y3) {};
    
    point2f p1_, p2_, p3_;
};

// transformation matrix (homogeneous 2d transformation matrix, 3 x 3)
class trans_matrix_t
{
public:
    
    trans_matrix_t(void)
    {
        m_.assign(3, vector<double>(3, 0));
        
        m_[0][0] = m_[1][1] = m_[2][2] = 1;
    };
    
    trans_matrix_t(const trans_matrix_t& mat)
    {
        m_.assign(mat.m_.begin(), mat.m_.end());
    };
    
    trans_matrix_t(const triangle_t& tri)
    {
        m_.assign(3, vector<double>(3, 0));
        
        m_[0][0] = tri.p1_.x_; m_[0][1] = tri.p2_.x_; m_[0][2] = tri.p3_.x_;
        
        m_[1][0] = tri.p1_.y_; m_[1][1] = tri.p2_.y_; m_[1][2] = tri.p3_.y_;
        
        m_[2].assign(3, 1);
    };
    
    // inverting a matrix
    trans_matrix_t invert(void)
    {
        trans_matrix_t mat;
        
        double determinant =  m_[0][0] * (m_[1][1] * m_[2][2] - m_[1][2] * m_[2][1])
                            - m_[0][1] * (m_[1][0] * m_[2][2] - m_[1][2] * m_[2][0])
                            + m_[0][2] * (m_[1][0] * m_[2][1] -m_[1][1] * m_[2][0]);
        
        if (0 == determinant)
        {
            cout << "determinant zero !" << endl; return trans_matrix_t();
        }
        
        mat[0][0] = (m_[1][1] * m_[2][2] - m_[1][2] * m_[2][1]) / determinant;
        mat[0][1] = (m_[0][2] * m_[2][1] - m_[0][1] * m_[2][2]) / determinant;
        mat[0][2] = (m_[0][1] * m_[1][2] - m_[0][2] * m_[1][1]) / determinant;
        
        mat[1][0] = (m_[1][2] * m_[2][0] - m_[1][0] * m_[2][2]) / determinant;
        mat[1][1] = (m_[0][0] * m_[2][2] - m_[0][2] * m_[2][0]) / determinant;
        mat[1][2] = (m_[0][2] * m_[1][0] - m_[0][0] * m_[1][2]) / determinant;
        
        mat[2][0] = (m_[1][0] * m_[2][1] - m_[1][1] * m_[2][0]) / determinant;
        mat[2][1] = (m_[0][1] * m_[2][0] - m_[0][0] * m_[2][1]) / determinant;
        mat[2][2] = (m_[0][0] * m_[1][1] - m_[0][1] * m_[1][0]) / determinant;
        
        return mat;
    };
    
    // get transpose of a matrix
    trans_matrix_t tranpose(void)
    {
        trans_matrix_t mat;
        
        for (int i = 0; i < 3; i ++)
            for (int k = 0; k < 3; k ++)
                mat[i][k] = m_[k][i];
        
        return mat;
    };
    
    // overloaded [] to enable [][]
    vector<double>& operator [](const int idx)
    {
        return m_[idx];
    };
    
    // matrix multiplication
    trans_matrix_t operator *(trans_matrix_t& mat)
    {
        trans_matrix_t ret_mat;
        
        for (int i = 0; i < 3; i ++)
        {
            for (int k = 0; k < 3; k ++)
            {
                ret_mat[i][k] = 0;
                
                for (int t = 0; t < 3; t ++)
                    ret_mat[i][k] += m_[i][t] * mat[t][k];
            }
        }
        
        return ret_mat;
    };
    
    vector<vector<double>> m_;
};

// triangle object (triangle data, transformation matrix, color)
class tri_object
{
public:
    
    tri_object(void) {};
    
    tri_object(const triangle_t& tri) : tri_(tri) {};
    
    tri_object(const triangle_t& tri,
               const trans_matrix_t& trans_mat)
                : tri_(tri), trans_mat_(trans_mat) {};
    
    tri_object(const triangle_t& tri,
               const trans_matrix_t& trans_mat,
               const color3f& color)
                : tri_(tri), trans_mat_(trans_mat), color_(color) {};
    
    triangle_t tri_;
    trans_matrix_t trans_mat_;
    color3f color_;
};

typedef complex<double> complex_t;          // position in complex plane
typedef vector<complex_t> complexlist_t;    // array of complex numbers
typedef vector<uint8_t> image_t;            // texture pixel data (unsigned char)
typedef vector<color3f> colorlist_t;        // color array
typedef vector<tri_object> objectlist_t;    // rendered object type
typedef vector<point2f> pointlist_t;        // point list type
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// global controlling variables

enumAssignmentType ASSIGNMENT_TYPE = MANDABROT_SET;     // current view type

colorlist_t COLOR_LIST;                 // global color list

int TEXTURE_WIDTH = 640;                // Mandabrot/Julia set picture size (in pixels)
int TEXTURE_HEIGHT = 480;

complex_t BOTTOM_LEFT_COMPLEX(-2.25, -1.5);  // window size of complex plane
complex_t TOP_RIGHT_COMPLEX(0.75, 1.5);

complex_t MANDABROT_BOTTOM_LEFT(-2.25, -1.5);  // window size of complex plane
complex_t MANDABROT_TOP_RIGHT(0.75, 1.5);         // for mandabrot set

complex_t JULIA_BOTTOM_LEFT(-2.0, -2.0);  // window size of complex plane
complex_t JULIA_TOP_RIGHT(2.0, 2.0);         // for mandabrot set

complex_t JULIA_SET_C(-0.8, -0.156);           // julia set parameter c

double MAG_RATIO = 1.0;                     // magnification ratio

image_t DISPLAY_IMAGE;                      // Mandabrot/Julia set image data

objectlist_t OBJECT_LIST;                   // list of objects to be drawn (triangles)

int IFS_RECURSION_DEPTH = 4;                // IFS recursion depth (k)

const double _NAN_ = -10000;
point2f OLD_MOUSE_POS(_NAN_, _NAN_);        // last mouse position

int WINDOW_WIDTH = 640;                     // window size
int WINDOW_HEIGHT = 480;

pointlist_t POINT_LIST_BUFFER;              // buffer for storing clicked vertices
//-------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------
void zoom(const bool is_zoom_in, complex_t& bottom_left, complex_t& top_right)
{
    double delta_x = top_right.real() - bottom_left.real();
    double delta_y = top_right.imag() - bottom_left.imag();
    
    double ratio = is_zoom_in ? 0.5 : -1.0;
    
    bottom_left += complex_t(delta_x * ratio / 2, delta_y * ratio / 2);
    top_right -= complex_t(delta_x * ratio / 2, delta_y * ratio / 2);
}

//
void zoom(const bool is_zoom_in, const enumAssignmentType& assignment_type)
{
    complex_t bottom_left;
    complex_t top_right;
    
    if (MANDABROT_SET == assignment_type)
    {
        bottom_left = MANDABROT_BOTTOM_LEFT;
        top_right = MANDABROT_TOP_RIGHT;
        
        zoom(is_zoom_in, bottom_left, top_right);
        
        MANDABROT_BOTTOM_LEFT = bottom_left;
        MANDABROT_TOP_RIGHT = top_right;
    }
    else if (JULIA_SET == assignment_type)
    {
        bottom_left = JULIA_BOTTOM_LEFT;
        top_right = JULIA_TOP_RIGHT;
        
        zoom(is_zoom_in, bottom_left, top_right);
        
        JULIA_BOTTOM_LEFT = bottom_left;
        JULIA_TOP_RIGHT = top_right;
    }
};

void move(const complex_t& vect, complex_t& bottom_left, complex_t& top_right)
{
    bottom_left += vect;
    top_right += vect;
}

//
void move(const complex_t& vect, const enumAssignmentType& assignment_type)
{
    complex_t bottom_left;
    complex_t top_right;
    
    if (MANDABROT_SET == assignment_type)
    {
        bottom_left = MANDABROT_BOTTOM_LEFT;
        top_right = MANDABROT_TOP_RIGHT;
        
        move(vect, bottom_left, top_right);
        
        MANDABROT_BOTTOM_LEFT = bottom_left;
        MANDABROT_TOP_RIGHT = top_right;
    }
    else if (JULIA_SET == assignment_type)
    {
        bottom_left = JULIA_BOTTOM_LEFT;
        top_right = JULIA_TOP_RIGHT;
        
        move(vect, bottom_left, top_right);
        
        JULIA_BOTTOM_LEFT = bottom_left;
        JULIA_TOP_RIGHT = top_right;
    }
};

//
void relocate(const double minx, const double maxx, const double miny, const double maxy,
              const enumAssignmentType& assignment_type)
{
    if (MANDABROT_SET == assignment_type)
    {
        MANDABROT_BOTTOM_LEFT = complex_t(minx, miny);
        MANDABROT_TOP_RIGHT = complex_t(maxx, maxy);
    }
    else if (JULIA_SET == assignment_type)
    {
        JULIA_BOTTOM_LEFT = complex_t(minx, miny);
        JULIA_TOP_RIGHT = complex_t(maxx, maxy);
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// TODO: Provide your solution within this area

//-------------------------------------------------------------------------------------
// TODO: use our provided color table or setup your own (optional)
//
// OUTPUT: color table represented as std::vector<color3f>
//
// HINT: check http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
//
void setupColorTable(colorlist_t& color_list)
{
    double color[] = {0.670588, 0, 0,
        0.682353, 0, 0,
        0.694118, 0, 0,
        0.705882, 0, 0,
        0.717647, 0, 0,
        0.729412, 0, 0,
        0.741176, 0, 0,
        0.752941, 0, 0,
        0.764706, 0, 0,
        0.776471, 0, 0,
        0.788235, 0, 0,
        0.8, 0, 0,
        0.811765, 0, 0,
        0.823529, 0, 0,
        0.835294, 0, 0,
        0.847059, 0, 0,
        0.858824, 0, 0,
        0.870588, 0, 0,
        0.882353, 0, 0,
        0.894118, 0, 0,
        0.905882, 0, 0,
        0.917647, 0, 0,
        0.929412, 0, 0,
        0.937255, 0, 0,
        0.94902, 0, 0,
        0.960784, 0, 0,
        0.972549, 0, 0,
        0.984314, 0, 0,
        0.996078, 0, 0,
        1, 0.00784314, 0,
        1, 0.0196078, 0,
        1, 0.0313725, 0,
        1, 0.0431373, 0,
        1, 0.054902, 0,
        1, 0.0666667, 0,
        1, 0.0784314, 0,
        1, 0.0901961, 0,
        1, 0.101961, 0,
        1, 0.113725, 0,
        1, 0.12549, 0,
        1, 0.137255, 0,
        1, 0.14902, 0,
        1, 0.160784, 0,
        1, 0.172549, 0,
        1, 0.184314, 0,
        1, 0.196078, 0,
        1, 0.207843, 0,
        1, 0.219608, 0,
        1, 0.231373, 0,
        1, 0.243137, 0,
        1, 0.254902, 0,
        1, 0.266667, 0,
        1, 0.278431, 0,
        1, 0.290196, 0,
        1, 0.301961, 0,
        1, 0.313725, 0,
        1, 0.32549, 0,
        1, 0.337255, 0,
        1, 0.34902, 0,
        1, 0.360784, 0,
        1, 0.372549, 0,
        1, 0.384314, 0,
        1, 0.396078, 0,
        1, 0.407843, 0,
        1, 0.419608, 0,
        1, 0.431373, 0,
        1, 0.443137, 0,
        1, 0.454902, 0,
        1, 0.466667, 0,
        1, 0.478431, 0,
        1, 0.490196, 0,
        1, 0.501961, 0,
        1, 0.513725, 0,
        1, 0.52549, 0,
        1, 0.537255, 0,
        1, 0.54902, 0,
        1, 0.560784, 0,
        1, 0.572549, 0,
        1, 0.584314, 0,
        1, 0.596078, 0,
        1, 0.607843, 0,
        1, 0.619608, 0,
        1, 0.631373, 0,
        1, 0.643137, 0,
        1, 0.654902, 0,
        1, 0.666667, 0,
        1, 0.678431, 0,
        1, 0.690196, 0,
        1, 0.701961, 0,
        1, 0.713725, 0,
        1, 0.72549, 0,
        1, 0.737255, 0,
        1, 0.74902, 0,
        1, 0.760784, 0,
        1, 0.772549, 0,
        1, 0.784314, 0,
        1, 0.796078, 0,
        1, 0.807843, 0,
        1, 0.819608, 0,
        1, 0.831373, 0,
        1, 0.843137, 0,
        1, 0.854902, 0,
        1, 0.866667, 0,
        1, 0.878431, 0,
        1, 0.890196, 0,
        1, 0.901961, 0,
        1, 0.913725, 0,
        1, 0.92549, 0,
        1, 0.937255, 0,
        1, 0.94902, 0,
        1, 0.960784, 0,
        1, 0.972549, 0,
        1, 0.984314, 0,
        1, 0.996078, 0,
        0.992157, 1, 0,
        0.968627, 1, 0.00392157,
        0.945098, 1, 0.00784314,
        0.921569, 1, 0.0117647,
        0.901961, 1, 0.0156863,
        0.878431, 1, 0.0196078,
        0.854902, 1, 0.0235294,
        0.831373, 1, 0.027451,
        0.807843, 1, 0.0313725,
        0.788235, 1, 0.0352941,
        0.764706, 1, 0.0392157,
        0.741176, 1, 0.0431373,
        0.717647, 1, 0.0470588,
        0.694118, 1, 0.054902,
        0.67451, 1, 0.0588235,
        0.65098, 1, 0.0627451,
        0.627451, 1, 0.0666667,
        0.603922, 1, 0.0705882,
        0.580392, 1, 0.0745098,
        0.556863, 1, 0.0784314,
        0.537255, 1, 0.0823529,
        0.513725, 1, 0.0862745,
        0.490196, 1, 0.0901961,
        0.466667, 1, 0.0941176,
        0.443137, 1, 0.0980392,
        0.423529, 1, 0.101961,
        0.4, 1, 0.105882,
        0.376471, 1, 0.109804,
        0.352941, 1, 0.113725,
        0.341176, 0.996078, 0.12549,
        0.329412, 0.996078, 0.137255,
        0.313725, 0.992157, 0.14902,
        0.301961, 0.988235, 0.156863,
        0.290196, 0.988235, 0.168627,
        0.278431, 0.984314, 0.180392,
        0.266667, 0.984314, 0.192157,
        0.254902, 0.980392, 0.203922,
        0.239216, 0.980392, 0.211765,
        0.227451, 0.976471, 0.223529,
        0.215686, 0.976471, 0.235294,
        0.203922, 0.972549, 0.247059,
        0.192157, 0.972549, 0.258824,
        0.180392, 0.968627, 0.270588,
        0.164706, 0.968627, 0.278431,
        0.152941, 0.964706, 0.290196,
        0.141176, 0.964706, 0.301961,
        0.129412, 0.960784, 0.313725,
        0.117647, 0.960784, 0.32549,
        0.105882, 0.956863, 0.333333,
        0.0901961, 0.952941, 0.345098,
        0.0784314, 0.952941, 0.356863,
        0.0666667, 0.94902, 0.368627,
        0.054902, 0.94902, 0.380392,
        0.0431373, 0.945098, 0.392157,
        0.0313725, 0.945098, 0.4,
        0.0156863, 0.941176, 0.411765,
        0.00392157, 0.941176, 0.423529,
        0, 0.929412, 0.439216,
        0, 0.905882, 0.458824,
        0, 0.886275, 0.478431,
        0, 0.862745, 0.498039,
        0, 0.839216, 0.517647,
        0, 0.819608, 0.537255,
        0, 0.796078, 0.556863,
        0, 0.772549, 0.580392,
        0, 0.752941, 0.6,
        0, 0.729412, 0.619608,
        0, 0.709804, 0.639216,
        0, 0.686275, 0.658824,
        0, 0.662745, 0.678431,
        0, 0.643137, 0.698039,
        0, 0.619608, 0.717647,
        0, 0.596078, 0.737255,
        0, 0.576471, 0.756863,
        0, 0.552941, 0.780392,
        0, 0.533333, 0.8,
        0, 0.509804, 0.819608,
        0, 0.486275, 0.839216,
        0, 0.466667, 0.858824,
        0, 0.443137, 0.878431,
        0, 0.419608, 0.898039,
        0, 0.4, 0.917647,
        0, 0.376471, 0.937255,
        0, 0.356863, 0.960784,
        0, 0.333333, 0.980392,
        0, 0.309804, 0.996078,
        0, 0.301961, 0.992157,
        0, 0.290196, 0.984314,
        0, 0.278431, 0.976471,
        0, 0.266667, 0.968627,
        0, 0.254902, 0.964706,
        0, 0.247059, 0.956863,
        0, 0.235294, 0.94902,
        0, 0.223529, 0.941176,
        0, 0.211765, 0.937255,
        0, 0.2, 0.929412,
        0, 0.192157, 0.921569,
        0, 0.180392, 0.913725,
        0, 0.168627, 0.909804,
        0, 0.156863, 0.901961,
        0, 0.145098, 0.894118,
        0, 0.137255, 0.886275,
        0, 0.12549, 0.882353,
        0, 0.113725, 0.87451,
        0, 0.101961, 0.866667,
        0, 0.0901961, 0.858824,
        0, 0.0823529, 0.854902,
        0, 0.0705882, 0.847059,
        0, 0.0588235, 0.839216,
        0, 0.0470588, 0.831373,
        0, 0.0352941, 0.827451,
        0, 0.027451, 0.819608,
        0, 0.0156863, 0.811765,
        0, 0.00392157, 0.803922};
    
    color_list.assign(228, color3f());
    for (size_t i = 0; i < color_list.size(); i ++)
    {
        color_list[i] = color3f(color[3*i], color[3*i+1], color[3*i+2]);
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// TODO: calculate Mandabrot set
//
// INPUT: texture image size, window area on complex plane, color table
// OUTPUT: Mandabrot set image within selected window area
//         image type std::vector<uint8_t>(image_width * image_height * 4)
//         pixels stored in row major order, pixel data: (red, green, blue, 255)
//
//
// image_width, image_height: size of the mandabrot set image
//
// c_bottom_left, c_top_right: bottom-left, top-right corner of the window
//
// color_list: color table
//
// image_t& image
//
void setMandabrotSet(const int image_width, const int image_height,
                     const complex_t& c_bottom_left, const complex_t& c_top_right,
                     const colorlist_t& color_list, image_t& image)
{
    // TODO: replace the example with your solution here
    
    /*size_t color_count = color_list.size();
    
    int image_size = image_width * image_height;

    image.assign(image_size * 4, 0);
    for (int x = 0; x < image_width; x ++)
    {
        for (int y = 0; y < image_height; y ++)
        {
            int pixel_idx = y * image_width + x;
            
            color3f color = color_list[(x + y) % color_count];
            
            image[pixel_idx * 4] = (uint8_t)(color.r_ * 255);
            image[pixel_idx * 4 + 1] = (uint8_t)(color.g_ * 255);
            image[pixel_idx * 4 + 2] = (uint8_t)(color.b_ * 255);
            image[pixel_idx * 4 + 3] = 255;
        }
    }*/
    
    int image_size = image_width * image_height;
    size_t color_count = color_list.size();
    image.assign(image_size * 4, 0);
    double da = (c_top_right.real() - c_bottom_left.real()) / image_width;
    double db = (c_top_right.imag()- c_bottom_left.imag()) / image_height;
    double x = 0.0; //real part of z
    double y = 0.0; //complex part of z
    double a = 0.0; //real part of c
    double b = 0.0; // complex part of c
    int n = 1;
    
    for(int p = 0; p < image_width; p++){
        a = c_bottom_left.real() + p * da;
        for(int q = 0; q < image_height;q++){
           b = c_bottom_left.imag() + q * db;
           y = b;
           x = a;
           n = 0;
           int pixel_idx = q * image_width + p;
            color3f color = color_list[n];

            while ((x*x + y*y) < 8 && n < color_count) {
                double nextx = x * x - y * y + a;
                double nexty = 2 * x * y + b;
                x = nextx;
                y = nexty;
                n++;
            }
            
            if( n == color_count && (x*x + y*y) < 8){
                //color pixel (p,q) with color 0, point in the set
                image[pixel_idx * 4] = (uint8_t)(color.r_ * 0);
                image[pixel_idx * 4 + 1] = (uint8_t)(color.g_ * 0);
                image[pixel_idx * 4 + 2] = (uint8_t)(color.b_ * 0);
                image[pixel_idx * 4 + 3] = 255;
            }else{
                if(n == 0){
                    color = color_list[n];
                }else{
                    color = color_list[n-1];
                }
                   //color pixel(p,q) with color n;
                    image[pixel_idx * 4] = (uint8_t)(color.r_ * 255);
                    image[pixel_idx * 4 + 1] = (uint8_t)(color.g_ * 255);
                    image[pixel_idx * 4 + 2] = (uint8_t)(color.b_ * 255);
                    image[pixel_idx * 4 + 3] = 255;
            }
        }
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// TODO: calculate Julia set
//
// INPUT: texture image size, window area on complex plane, c, color table
//
// OUTPUT: Julia set image within selected window area
//         image type std::vector<uint8_t>(image_width * image_height * 4)
//         pixels stored in row major order, pixel data: (red, green, blue, 255)
//
// image_width, image_height: size of the julia set image
//
// c_bottom_left, c_top_right: bottom-left, top-right corner of the window
//
// c: the controlling arguments c in Julia sets
//
// color_list: color table
//
// image_t& image
//
void setJuliaSet(const int image_width, const int image_height,
                 const complex_t& c_bottom_left, const complex_t& c_top_right,
                 const complex_t& c, const colorlist_t& color_list, image_t& image)
{
    // TODO: replace the example with your solution here
    
    /*size_t color_count = color_list.size();
    
    int image_size = image_width * image_height;
    
    image.assign(image_size * 4, 0);
    for (int x = 0; x < image_width; x ++)
    {
        for (int y = 0; y < image_height; y ++)
        {
            int pixel_idx = y * image_width + x;
            
            color3f color = color_list[(x - y) % color_count];
            
            image[pixel_idx * 4] = (uint8_t)(color.r_ * 255);
            image[pixel_idx * 4 + 1] = (uint8_t)(color.g_ * 255);
            image[pixel_idx * 4 + 2] = (uint8_t)(color.b_ * 255);
            image[pixel_idx * 4 + 3] = 255;
        }
    }*/
   // complex_t JULIA_BOTTOM_LEFT(-2.0, -2.0);  // window size of complex plane
   // complex_t JULIA_TOP_RIGHT(2.0, 2.0);         // for mandabrot set
    int image_size = image_width * image_height;
    size_t color_count = color_list.size();
    image.assign(image_size * 4, 0);
    double dx = (c_top_right.real() - c_bottom_left.real()) / image_width;
    double dy = (c_top_right.imag()- c_bottom_left.imag()) / image_height;
    double x = 0.0; //real part of z
    double y = 0.0; //complex part of z
    double a = c.real();//real part of c
    double b = c.imag(); // complex part of c
    double initx;
    double inity;
    int n = 1;
    
    for(int p = 0; p < image_width; p++){
        initx = c_bottom_left.real() + p * dx;
        for(int q = 0; q < image_height;q++){
            inity = c_bottom_left.imag() + q * dy;
            y = inity;
            x = initx;
            n = 0;
            int pixel_idx = q * image_width + p;
            color3f color = color_list[n];
            
            while ((x*x + y*y) < 4 && n < color_count) {
                double nextx = x * x - y * y + a;
                double nexty = 2 * x * y + b;
                x = nextx;
                y = nexty;
                n++;
            }
            
            if( n == color_count && (x*x + y*y) < 4){
                //color pixel (p,q) with color 0, point in the set
                image[pixel_idx * 4] = (uint8_t)(color.r_ * 0);
                image[pixel_idx * 4 + 1] = (uint8_t)(color.g_ * 0);
                image[pixel_idx * 4 + 2] = (uint8_t)(color.b_ * 0);
                image[pixel_idx * 4 + 3] = 255;
            }else{
                //color pixel(p,q) with color n;
                if(n == 0){
                    color = color_list[n];
                }else{
                    color = color_list[n-1];
                }
                image[pixel_idx * 4] = (uint8_t)(color.r_ * 255);
                image[pixel_idx * 4 + 1] = (uint8_t)(color.g_ * 255);
                image[pixel_idx * 4 + 2] = (uint8_t)(color.b_ * 255);
                image[pixel_idx * 4 + 3] = 255;
            }
        }
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// TODO: calculate the transformation matrix from triangle tri1 to triangle tri2
trans_matrix_t calcTransformationMatrix(const triangle_t& tri1, const triangle_t& tri2)
{
    trans_matrix_t mat1 = trans_matrix_t(tri1).invert();
    
    trans_matrix_t mat2(tri2);
    
    trans_matrix_t trans_mat = mat2 * mat1;
    
    return trans_mat;
};
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// add triangles when the size of input points > 3
// x_screen, y_screen: screen coordinates
// point_list: stored all clicked points
void pointClicked(const double x_screen, const double y_screen,
                  const pointlist_t& point_list)
{
    POINT_LIST_BUFFER.push_back(point2f(x_screen, y_screen));
    
    if (0 == POINT_LIST_BUFFER.size() % 3)
    {
        triangle_t tri_ref(POINT_LIST_BUFFER[0], POINT_LIST_BUFFER[1], POINT_LIST_BUFFER[2]);
        
        size_t size = POINT_LIST_BUFFER.size();
        
        triangle_t tri(POINT_LIST_BUFFER[size-3], POINT_LIST_BUFFER[size-2], POINT_LIST_BUFFER[size-1]);
        
        trans_matrix_t mat = calcTransformationMatrix(tri_ref, tri);
        
        tri_object obj(tri, mat);
        
        OBJECT_LIST.push_back(obj);
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// draw a triangle using OpenGL
void drawTriangle(const tri_object& tri_obj)
{
    triangle_t tri = tri_obj.tri_;

    color3f color = tri_obj.color_;
    glDisable(GL_TEXTURE_2D);
    
    glColor3f(color.r_, color.g_, color.b_);
    
    glBegin(GL_LINE_LOOP);

    
    glVertex3f(tri.p1_.x_, tri.p1_.y_, 0);
    glVertex3f(tri.p2_.x_, tri.p2_.y_, 0);
    glVertex3f(tri.p3_.x_, tri.p3_.y_, 0);
    
    glEnd();
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// TODO: recursively draw triangles using IFS
//
// INPUT: triangle list and recursion depth k
//        each element of triangle list stores a triangle data,
//        2D homogeneous transformation matrix and color
//
// OUTPUT: recursively draw triangles using OpenGL
//
int countcolor = 10;
void drawIfsTriangles(objectlist_t& triangle_list, const int k,int flag)
{
    
    // TODO: replace the following example with your solution
    // NOTE: the following code is only for demonstration,
    //          not related to the solution, delete them after you read them
    if(flag == 1){
    triangle_t tri1(point2f(-1.0 / 4, 0), point2f(1.0 / 4, 0), point2f(0, 2.0 / 4));
    triangle_t tri2(point2f(-0.5 / 4, 1.0 / 4), point2f(0.5 / 4, 1.0 / 4), point2f(0, 2.0 / 4));
    triangle_t tri3(point2f(-1.0 / 4, 0), point2f(0, 0), point2f(-0.5 / 4, 1.0 / 4));
    triangle_t tri4(point2f(0, 0), point2f(1.0 / 4, 0), point2f(0.5 / 4, 1.0 / 4));
    
    tri_object tri_obj1(tri1);
    tri_object tri_obj2(tri2);
    tri_object tri_obj3(tri3);
    tri_object tri_obj4(tri4);
    
    tri_obj1.trans_mat_ = calcTransformationMatrix(tri1, tri1);
    tri_obj2.trans_mat_ = calcTransformationMatrix(tri1, tri2);
    tri_obj3.trans_mat_ = calcTransformationMatrix(tri1, tri3);
    tri_obj4.trans_mat_ = calcTransformationMatrix(tri1, tri4);
    
    tri_obj1.color_ = COLOR_LIST[50];
    tri_obj2.color_ = COLOR_LIST[100];
    tri_obj3.color_ = COLOR_LIST[160];
    tri_obj4.color_ = COLOR_LIST[200];
    
    OBJECT_LIST.push_back(tri_obj1);
    OBJECT_LIST.push_back(tri_obj2);
    OBJECT_LIST.push_back(tri_obj3);
    OBJECT_LIST.push_back(tri_obj4);
    }
    
    /*for (int i = 0; i < OBJECT_LIST.size(); i ++)
        drawTriangle(OBJECT_LIST[i]);*/
    
    // TODO: insert your solution here, recursively draw triangles
    glMatrixMode(GL_MODELVIEW);
    //glScaled(0.2, 0.2, 0.2);
    int times = (unsigned int)triangle_list.size();
    if(k > 0){
        for(int count = 0; count < times;count++){
            glPushMatrix();
            triangle_list[count].trans_mat_ = calcTransformationMatrix(triangle_list[0].tri_, triangle_list[count].tri_);
            double *arr = new double [16];
            trans_matrix_t t = triangle_list[count].trans_mat_;
                    arr[0] = t[0][0];
                    arr[1] = t[1][0];
                    arr[2] = 0;
                    arr[3] = t[2][0];
            
                    arr[4] = t[0][1];
                    arr[5] = t[1][1];
                    arr[6] = 0;
                    arr[7] = t[2][1];
            
                    arr[8] = 0;
                    arr[9] = 0;
                    arr[10] = 0;
                    arr[11] = 0;
            
                    arr[12] = t[0][2];
                    arr[13] = t[1][2];
                    arr[14] = 0;
                    arr[15] = t[2][2];
            glMultMatrixd(arr);//current * arr
            drawIfsTriangles(triangle_list,k-1,0);
            glPopMatrix();
        }
    }else{
        /*GLdouble mat[16];
        glGetDoublev(GL_MODELVIEW_MATRIX, mat);
        trans_matrix_t m;

        m[0][0] = mat[0];
        m[1][0] = mat[1];
        m[2][0] = mat[3];
        m[0][1] = mat[4];
        m[1][1] = mat[5];
        m[2][1] = mat[7];
        m[0][2] = mat[12];
        m[1][2] = mat[13];
        m[2][2] = mat[15];*/
        countcolor = 10 + countcolor;
        trans_matrix_t tri(triangle_list[0].tri_);
       /* trans_matrix_t tri =  m1;*/
    triangle_t t(point2f(tri[0][0],tri[1][0]), point2f(tri[0][1],tri[1][1]), point2f(tri[0][2],tri[1][2]));
        triangle_list[0].tri_ = t;
        tri_object tri_obj(t);
        if(countcolor > 200){ countcolor = 10;}
        tri_obj.color_ = COLOR_LIST[countcolor];
        drawTriangle(tri_obj);
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// update Mandabrot/Julia set
void updateData(void)
{
    if (MANDABROT_SET == ASSIGNMENT_TYPE)
    {
        BOTTOM_LEFT_COMPLEX = MANDABROT_BOTTOM_LEFT;
        TOP_RIGHT_COMPLEX = MANDABROT_TOP_RIGHT;
        
        setMandabrotSet(TEXTURE_WIDTH, TEXTURE_HEIGHT,
                        BOTTOM_LEFT_COMPLEX, TOP_RIGHT_COMPLEX,
                        COLOR_LIST, DISPLAY_IMAGE);
    }
    else if (JULIA_SET == ASSIGNMENT_TYPE)
    {
        BOTTOM_LEFT_COMPLEX = JULIA_BOTTOM_LEFT;
        TOP_RIGHT_COMPLEX = JULIA_TOP_RIGHT;
        
        setJuliaSet(TEXTURE_WIDTH, TEXTURE_HEIGHT,
                    BOTTOM_LEFT_COMPLEX, TOP_RIGHT_COMPLEX, JULIA_SET_C,
                    COLOR_LIST, DISPLAY_IMAGE);
    }
    else if (IFS_TRIANGLES == ASSIGNMENT_TYPE)
    {
        ;
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// draw textures
void drawTextures(void)
{
    // current drawing color
    glColor3f(1.0, 1.0, 1.0);
    
    // enable 2D texture feature
    glEnable(GL_TEXTURE_2D);
    
    // texture ID
    GLuint tex_id = 0;
    
    // generate texture ID
    glGenTextures(1, &tex_id);
    
    // use texture represented by tex_id
    glBindTexture(GL_TEXTURE_2D, tex_id);
    
    // set texture parameters
    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    
    //
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    
    //
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
    
    // set texture data
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, TEXTURE_WIDTH, TEXTURE_HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, &DISPLAY_IMAGE[0]);
    
    // drawing, try different primitives
    glBegin(GL_QUADS);
    
    double square_size = 0.5;
    
    glTexCoord2f(0, 0); glVertex3f(-square_size, -square_size, 0.0);
    glTexCoord2f(1, 0); glVertex3f(square_size, -square_size, 0.0);
    glTexCoord2f(1, 1); glVertex3f(square_size, square_size, 0.0);
    glTexCoord2f(0, 1); glVertex3f(-square_size, square_size, 0.0);
    
    glEnd();
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// draw triangles
void drawTriangles(void)
{
    drawIfsTriangles(OBJECT_LIST, IFS_RECURSION_DEPTH,1);
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// update scene
void updateRendering()
{
    if (MANDABROT_SET == ASSIGNMENT_TYPE || JULIA_SET == ASSIGNMENT_TYPE)
        drawTextures();
    else if (IFS_TRIANGLES == ASSIGNMENT_TYPE)
        drawTriangles();
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// resize callback function (called each time when the window is resized)
// width / height: current window width / height passed by GLUT
void reshape(int width, int height)
{
    WINDOW_WIDTH = width; WINDOW_HEIGHT = height;
    
    updateData();
    
    //cout << "resized as new window size: ( " << width << ", " << height << " )" << endl;
    
    // this call just remind us what the current matrix is
    // delete this line when you are confident about this
    glMatrixMode(GL_MODELVIEW);
    
    // reset current view port to be the same size as the resized window
    glViewport(0, 0, (GLsizei)width, (GLsizei)height);
    
    // switch to projection matrix
    glMatrixMode(GL_PROJECTION);
    
    // reset the projection matrix to be identity matrix
    glLoadIdentity();
    
    // FOV angles (in degrees, 30 left, 30 right), width / height, near & far plane
    gluPerspective(60, (GLfloat)width / (GLfloat)height, 1.0, 100.0);
    
    // DO NOT FORGET THIS LINE, changing current matrix to be the Model-View Matrix
    // (leave the mode of Projection Matrix)
    glMatrixMode(GL_MODELVIEW);
    
    glutPostRedisplay();
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// display callback function (called each time during the rendering loop)
void display(void)
{
    //cout << "redisplaying" << endl;
    
    // realize what matrix is currently using
    // delete this line when you are confident about this
    glMatrixMode(GL_MODELVIEW);
    
    // set clearing color to be (0.5, 0.5, 0.5, 1.0) - grey
    glClearColor(0.5, 0.5, 0.5, 1.0);
    
    // reset the matrix
    glLoadIdentity();
    
    // clear the buffer using present clearing color
    glClear(GL_COLOR_BUFFER_BIT);
    
    // translate the scene along direction (0, 0, -5.0) - pushing them forward
    // a trick to move the camera backward
    glTranslatef(0.0f, 0.0f, -1.0f);
    
    // draw our things here
    updateRendering();
    
    // flush the painted things
    glFlush();
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// calculate world coordinate from screen position
void getWorldCoord(const int x_screen, const int y_screen, double& x, double& y)
{
    //
    
    GLdouble modelview[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    
    GLdouble projection[16];
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    
    GLfloat winZ;
    glReadPixels(x_screen, y_screen, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);

    double z;
    gluUnProject(x_screen, y_screen, winZ, modelview, projection, viewport, &x, &y, &z);
    
    
    
//    GLfloat winX = 0.0, winY = 0.0, winZ = 0.0;
//    
//    
//    GLdouble posX = 0.0, posY = 0.0, posZ = 0.0;
//    winX = (float)x;
//    winY = (float)OGLMviewport[3] - (float)y;   // invert winY so that down lowers value
//    glReadPixels( x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
//    gluUnProject( winX, winY, winZ, OGLMmodelview, OGLMprojection, OGLMviewport, &posX, &posY, &posZ);

    
    cout << "returned world position: ( " << x << ", " << y << " )" << endl << endl;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// mouse pressing callback function (called when certain mouse event happens)
// key: character representing that key
// x, y: the mouse position when the key is pressed
void mousePress(int button, int state, int x, int y)
{
    //cout << "mouse position: (" << x << ", " << y << ")" << endl;

    //cout << "mouse clicked" << endl;
    
    int y_temp = WINDOW_HEIGHT - 1 - y;
    
    if (MANDABROT_SET == ASSIGNMENT_TYPE || JULIA_SET == ASSIGNMENT_TYPE)
    {
        if (GLUT_LEFT_BUTTON == button && GLUT_DOWN == state)
        {
            //if ( _NAN_ == OLD_MOUSE_POS.x_ && _NAN_ == OLD_MOUSE_POS.y_ )
                OLD_MOUSE_POS = point2f(x, y_temp);
        }
        else if (GLUT_LEFT_BUTTON == button && GLUT_UP == state)
        {
            if (_NAN_ != OLD_MOUSE_POS.x_ && _NAN_ != OLD_MOUSE_POS.y_)
            {
                double delta_x = (x - OLD_MOUSE_POS.x_);
                double delta_y = (y_temp - OLD_MOUSE_POS.y_);
                
                double mag = sqrt(delta_x * delta_x + delta_y * delta_y);
                
                if(mag == 0) return;
                delta_x /= mag;
                delta_y /= mag;
                
                move(complex_t(-0.2 * delta_x, -0.2 * delta_y) / MAG_RATIO, ASSIGNMENT_TYPE);
                
                updateData();
                
                glutPostRedisplay();
                
                OLD_MOUSE_POS = point2f(_NAN_, _NAN_);
            }
        }
        
        if (JULIA_SET == ASSIGNMENT_TYPE
            && GLUT_MIDDLE_BUTTON == button && GLUT_DOWN == state)
        {
//            double x_coord, y_coord;
//            
//            // calculate world coordinate from mouse clicked points
//            getWorldCoord(x, y_temp, x_coord, y_coord);
//            
//            JULIA_SET_C = complex_t(x_coord, y_coord);
//            
//            updateData();
//            
//            glutPostRedisplay();
        }
    }
    else if (IFS_TRIANGLES == ASSIGNMENT_TYPE)
    {
        if (GLUT_LEFT_BUTTON == button && GLUT_DOWN == state)
        {
            double x_coord, y_coord;
            
            // calculate world coordinate from mouse clicked points
            getWorldCoord(x, y_temp, x_coord, y_coord);
            
            pointClicked(x_coord, y_coord, POINT_LIST_BUFFER);
            
            glutPostRedisplay();
        }
    }
}
//-------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------
// mouse moving callback function (called when mouse moves with button clicked)
// x, y: the mouse position
void mouseMove(int x, int y)
{
    //cout << "mouse moving: ( " << x << ", " << y << ")" << endl;
    
    // ...
}

//-------------------------------------------------------------------------------------
// key pressing callback function (called when a key is pressed)
// key: character representing that key
// x, y: the mouse position when the key is pressed
void keyPress(unsigned char key, int x, int y)
{
    //cout << "pressed key: \"" << key << "\", position: (" << x << ", " << y << ")" << endl;
    
    if (13 == key)
    {
        ASSIGNMENT_TYPE = static_cast<enumAssignmentType>((ASSIGNMENT_TYPE + 1) % 3);
        
        if (MANDABROT_SET == ASSIGNMENT_TYPE)
        {
            BOTTOM_LEFT_COMPLEX = MANDABROT_BOTTOM_LEFT;
            TOP_RIGHT_COMPLEX = MANDABROT_TOP_RIGHT;
            
            MAG_RATIO = 1.0;
        }
        else if (JULIA_SET == ASSIGNMENT_TYPE)
        {
            BOTTOM_LEFT_COMPLEX = JULIA_BOTTOM_LEFT;
            TOP_RIGHT_COMPLEX = JULIA_TOP_RIGHT;
            
            MAG_RATIO = 1.0;
        }
    }
    
    if (MANDABROT_SET == ASSIGNMENT_TYPE || JULIA_SET == ASSIGNMENT_TYPE)
    {
        if ('+' == key)
        {
            MAG_RATIO *= 1.5;
            zoom(true, ASSIGNMENT_TYPE);
        }
        else if ('-' == key)
        {
            MAG_RATIO /= 1.5;
            zoom(false, ASSIGNMENT_TYPE);
        }
    }
    else if (IFS_TRIANGLES == ASSIGNMENT_TYPE)
    {
        if ('+' == key)
        {
            IFS_RECURSION_DEPTH ++;
        }
        else if ('-' == key)
        {
            IFS_RECURSION_DEPTH --;
            
            if (IFS_RECURSION_DEPTH < 0) IFS_RECURSION_DEPTH = 0;
        }
    }
    
    updateData();
    
    glutPostRedisplay();
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// key releasing callback function (called when a key is released)
// key: character representing that key
// x, y: the mouse position when the key is pressed
void keyRelease(unsigned char key, int x, int y)
{
    //cout << "released key: \"" << key << "\", position: (" << x << ", " << y << ")" << endl;
    
    if ('w' == key)
    {
        // do things when 'w' is pressed
    }
    // ...
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// special key pressing callback function (called when a special key is pressed)
// specialKey: integer representing the special key
// x, y: the mouse position when the key is pressed
void keySpecialPress(int specialKey, int x, int y)
{
    //cout << "pressed special key: \"" << specialKey << "\", position: (" << x << ", " << y << ")" << endl;
    
    if (MANDABROT_SET == ASSIGNMENT_TYPE || JULIA_SET == ASSIGNMENT_TYPE)
    {
        if (GLUT_KEY_LEFT == specialKey)
        {
            move(complex_t(0.2, 0) / MAG_RATIO, ASSIGNMENT_TYPE);
        }
        else if (GLUT_KEY_RIGHT == specialKey)
        {
            move(complex_t(-0.2, 0) / MAG_RATIO, ASSIGNMENT_TYPE);
        }
        else if (GLUT_KEY_DOWN == specialKey)
        {
            move(complex_t(0, 0.2) / MAG_RATIO, ASSIGNMENT_TYPE);
        }
        else if (GLUT_KEY_UP == specialKey)
        {
            move(complex_t(0, -0.2) / MAG_RATIO, ASSIGNMENT_TYPE);
        }
        
        updateData();
        
        glutPostRedisplay();
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// special key releasing callback function (called when a special key is released)
// specialKey: integer representing the special key
// x, y: the mouse position when the key is pressed
void keySpecialRelease(int specialKey, int x, int y)
{
    //cout << "pressed special key: \"" << specialKey << "\", position: (" << x << ", " << y << ")" << endl;
    
    if (GLUT_KEY_LEFT == specialKey)
    {
        // do things when '<-' is pressed
    }
    // ...
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
// program entry point
int main(int argc, char ** argv)
{
    // setup color list
    setupColorTable(COLOR_LIST);
    
    // Initialize GLUT
    glutInit(&argc, argv);
    
    // set display mode (currently only single frame buffer)
    glutInitDisplayMode(GLUT_SINGLE);
    
    // window size & position
    glutInitWindowSize (WINDOW_WIDTH, WINDOW_HEIGHT);
    glutInitWindowPosition (100, 100);
    
    // create window
    glutCreateWindow ("Hello World - OpenGL");
    
    // register callback functions
    glutDisplayFunc(display);               // displaying
    glutReshapeFunc(reshape);               // resizing
    
    glutKeyboardFunc(keyPress);             // key pressing /releasing
    glutKeyboardUpFunc(keyRelease);
    
    glutMouseFunc(mousePress);            // mouse clicking
    glutMotionFunc(mouseMove);              // mouse moving
    
    glutSpecialFunc(keySpecialPress);       // special-key pressing / releasing
    glutSpecialUpFunc(keySpecialRelease);
    
    // looping
    glutMainLoop();
    
    return 0;
}
//-------------------------------------------------------------------------------------

