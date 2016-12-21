// #pragma once
#include <string.h>
#include <jni.h>
#include "opencv2/core/core.hpp"
#include <opencv2/imgproc/imgproc.hpp>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <android/log.h>
#include <map>

using namespace std;
using namespace cv;

/// Type Defs
// First int is ID, second int is counter and a contour
typedef std::map<int, std::pair<vector<Point>, int> > ContoursMap;

/// "Global" vars

RNG rng(12345);

int lastSecondTimeStamp = 0;

int secondsSinceStartUp = 0;

int outputTimeOutContours = 30;

int currentFrame = 0;

int numConsideredFrames = 7;	// Consider last x frames for calc new contours

int numberOfCannysPerFrame = 3;

bool bufferSwitch = true;

int remainingContourFrames = 0;

int remainingARImgFrames = 0;

int minSecsToShowARImg = 2;

int arImgFrames = 5;

int frameWidth;

int frameHeight;

int watermarkWidth;

int watermarkHeight;

float globalWaterMarkSize;

int globalWatermarkMode;

bool switchState = false;

/// Define structure dimensions and pos of zenith
// Feldherrnhalle
const float FELDHERRNHALLE_WIDTH = 106.4f;
const float FELDHERRNHALLE_HEIGHT = 46.0f;
const float FELDHERRNHALLE_ZENITH_PERC_POS_X = 0.5f;
const float FELDHERRNHALLE_ZENITH_PERC_POS_Y = 84.0f / 363.0f;
// Siegestor
const float SIEGESTOR_WIDTH = 125.0f;
const float SIEGESTOR_HEIGHT = 54.0f;
const float SIEGESTOR_PERC_POS_X = 0.5f;
const float SIEGESTOR_PERC_POS_Y = 265.0f / 546.0f;

// Init structure values
float archHeight = FELDHERRNHALLE_HEIGHT;
float archWidth = FELDHERRNHALLE_WIDTH;
float archZenithPercPosX = FELDHERRNHALLE_ZENITH_PERC_POS_X;
float archZenithPercPosY = FELDHERRNHALLE_ZENITH_PERC_POS_Y;
bool archHasBackground = true;

Mat originalRGBImg;
Mat arImage;
Mat resizedARImage;
Point touchedPoint;

int number_of_arch_widths = arImgFrames * 2;
int last_pushed_index = 0;

vector<int> last_arch_widths(number_of_arch_widths);

int old_middle_arch_width = 0;
vector<Point> old_middle_arch;
Rect old_boundingRec;

vector<vector<Point>> prevFoundArchContours;
vector<vector<vector<Point>>> switchBuffer(2);
vector<vector<Point>> searchContours;
ContoursMap globalContoursMap;

Point center;

// Define basic colors
Scalar white = Scalar(255, 255, 255);
Scalar black = Scalar(0, 0, 0);
Scalar red = Scalar(255, 0, 0);
Scalar green = Scalar(0, 255, 0);
Scalar blue = Scalar(0, 0, 255);

/// Method Declarations
void generateSearchContours();

void updateSettings(int width, int height);

void resetOldMiddleArchProps();

void applySobel(Mat& mGray, int kernel_size = 5);

void applyBlurAndCanny(Mat& mGray, double factor = 0.5);

void opticalDetectionDebug(Mat& mRgb, Mat& mGray);

void opticalDetection(Mat& mRgb, Mat& mGray);

bool similarWidth(vector<Point> contour1, vector<Point> contour2, float maxDeviation = 0.25f);

Rect createRect(Point p1, Point p2, int expandVal = 0);

Point getCenterOfContour(vector<Point> contour);

float getLongestContourSide(vector<Point> contour);

vector<Point> getMiddleArch(Mat mRgb, vector<vector<Point>> contours);

void mapARImageOnMiddleArch(Mat mRgb, vector<vector<Point>> contours, float middleArchDiameter);

void mapARImageOnZenithOfMiddleArch(Mat mRgb, vector<Point> middleArch, float middleArchDiameter);

Point getZenithOfContour(vector<Point> contour);

vector<Point> getAtPointedArch(Mat mRgb, vector<vector<Point>> contours);

bool contourIsArch(vector<Point> contour, bool drawPoints = false);

bool pointsAreAtOuterPos(Point point1, Point point2, vector<Point> contour);

bool pointIsAtBottom(Point point, vector<Point> contour);

float getSideLength(Point p1, Point p2);

float getOppositeAngleOfSideC(float sideA, float sideB, float sideC);

bool pointsAtSamePos(Point p1, Point p2, float maxDeviation = 0.2f);

bool contourCentersMatching(vector<Point> contour1, vector<Point> contour2, float maxBoundAreaDiff = 0.8f, float maxCenterDiff = 0.5f);

bool pointInRect(Point point, Rect rect);

vector<vector<Point>> eliminateDuplicateContours(vector<vector<Point>> contours);

vector<vector<Point>> extractContours(Mat grayImg, int minShapeArea = 50, int minSideLength = 20, int minVerts = 3, int maxVerts = 24);

vector<vector<Point>> pruneContours(Mat mRgb, vector<vector<Point>> contours, int minPoints = 4, float maxRotArea = 0.66f, double maxProportion = 4.0, bool onlySpecContours = true, double shapeAccuracy = 0.15);

vector<Point> extractArch(vector<Point> contour);

ContoursMap contoursToContoursMap(vector<vector<Point>> contours, int maxNumOfMisses = 3);

vector<vector<Point>> contoursMapToContours(ContoursMap contoursMap);

ContoursMap matchContoursMaps(vector<vector<Point>> newContours, ContoursMap oldContoursMap, float maxDeviation = 0.25f, float matchShapeDeviation = 0.05f, int maxNumOfMisses = 3);

vector<vector<Point>> matchContours(vector<vector<Point>> contours1, vector<vector<Point>> contours2, float maxDeviation = 0.25f, float matchShapeDeviation = 0.05f);

bool contoursRoughlyAtSamePos(vector<Point> contour1, vector<Point> contour2, float maxDeviation = 0.1f, int maxTotalDiffer = 1);

void setBlackPixelsTransparent(Mat image, Mat backgroundImg);

void changeAlphaValue(Mat image, Mat backgroundImg, float opaquePercentage = 0.4f);

void drawFoundContours(Mat srcDst, vector<vector<Point>> foundContours, int thickness = 2);

vector<Point> createArch(int radius, int numOfPointsInSemiCircle = 5, float scaleX = 1.0f, float scaleY = 1.0f, float circleScaleX = 0.33f, bool archExtension = true);

void drawBorder(Mat mRgb, int borderSize = 4);

/// Native Method Declarations
/*
To allow for overloading of functions, C++ uses something called ** name mangling **.
This means that function names are not the same in C++ as in plain C.
To inhibit this name mangling, you have to declare functions as extern “C”
*/
extern "C" {
	JNIEXPORT void JNICALL Java_com_tum_historicarguide_MainActivity_setupDetection(JNIEnv *env, jobject instance, jint width, jint height, jlong addrARImg);

	JNIEXPORT void JNICALL Java_com_tum_historicarguide_MainActivity_updateSettings(JNIEnv *env, jobject instance, jint width, jint height);

	JNIEXPORT void JNICALL Java_com_tum_historicarguide_MainActivity_setSwitchState(JNIEnv *env, jobject instance, jboolean switchState);

	JNIEXPORT void JNICALL Java_com_tum_historicarguide_MainActivity_nativeOpticalDetectionDebug(JNIEnv *env, jobject instance, jlong addrRgba, jfloat fps);

	JNIEXPORT void JNICALL Java_com_tum_historicarguide_MainActivity_nativeOpticalDetection(JNIEnv *env, jobject instance, jlong addrRgba, jfloat fps);

	JNIEXPORT void JNICALL Java_com_tum_historicarguide_MainActivity_nativeSetTouchPos(JNIEnv *env, jobject instance, jint xCoord, jint yCoord);
}// END extern "C"