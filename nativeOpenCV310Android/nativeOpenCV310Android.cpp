/**********************************
Java Native Interface library
**********************************/
#include "nativeOpenCV310Android.h"

extern "C" {
	JNIEXPORT void JNICALL Java_com_tum_historicarguide_MainActivity_setupDetection(JNIEnv *env, jobject instance, jint innerWidth, jint innerHeight, jint frameWidth, jint frameHeight, jlong addrARImg) {
		// Store images
		arImage = *(Mat*)addrARImg;
		
		std::ostringstream rowOstr;
		rowOstr << " ";
		__android_log_write(ANDROID_LOG_INFO, "setupDetection22222", rowOstr.str().c_str());

		center = Point((int)innerWidth / 2, (int)innerHeight / 2);
		// touchedPoint = center;

		// Init searchShapes
		searchContours.clear();
		generateSearchContours();

		updateSettings((int)innerWidth, (int)innerHeight, (int)frameWidth, (int)frameHeight);
	}

	// TODO: Not needed anymore
	/*
	JNIEXPORT void JNICALL Java_com_tum_historicarguide_MainActivity_updateSettings(JNIEnv *env, jobject instance, jint innerWidth, jint innerHeight) {
		updateSettings((int)innerWidth, (int)innerHeight);
	}
	*/

	JNIEXPORT void JNICALL Java_com_tum_historicarguide_MainActivity_setSwitchState(JNIEnv *env, jobject instance, jboolean switchState) {
		if (switchState) {
			archWidth = SIEGESTOR_WIDTH;
			archHeight = SIEGESTOR_HEIGHT;
			archZenithPercPosX = SIEGESTOR_PERC_POS_X;
			archZenithPercPosY = SIEGESTOR_PERC_POS_Y;
			archHasBackground = false;
		}
		else {
			archWidth = FELDHERRNHALLE_WIDTH;
			archHeight = FELDHERRNHALLE_HEIGHT;
			archZenithPercPosX = FELDHERRNHALLE_ZENITH_PERC_POS_X;
			archZenithPercPosY = FELDHERRNHALLE_ZENITH_PERC_POS_Y;
			archHasBackground = true;
		}

		resetOldMiddleArchProps();
	}

	Mat glob_mRgb;

	JNIEXPORT void JNICALL Java_com_tum_historicarguide_MainActivity_nativeOpticalDetectionDebug(JNIEnv *env, jobject instance, jlong innerRgba, jlong addrRgba, jfloat fps) {
		Mat& innerRgb = *(Mat*)innerRgba;
		glob_mRgb = *(Mat*)addrRgba;

		// cvtColor(mRgb, originalRGBImg, CV_RGBA2RGB);
		glob_mRgb.copyTo(originalRGBImg);

		Mat mGray;
		cvtColor(innerRgb, mGray, CV_RGBA2GRAY);

		Mat equalizedMat;
		equalizeHist(mGray, equalizedMat);

		arImgFrames = (int)fps;

		/*
		std::ostringstream rowOstr;
		rowOstr << innerRgb.cols;
		rowOstr << "/";
		rowOstr << innerRgb.rows;
		rowOstr << ", RGB SIZE: ";
		rowOstr << mRgb.cols;
		rowOstr << "/";j
		rowOstr << mRgb.rows;
		__android_log_write(ANDROID_LOG_INFO, "innerRgb SIZE: ", rowOstr.str().c_str());
		*/

		opticalDetectionDebug(innerRgb, equalizedMat);

		// glob_mRgb = originalRGBImg;

		// TODO: Not needed anymore
		/*		// Draw touch
		if (touchedPoint.x != 0 || touchedPoint.y != 0) {
			circle(innerRgb, touchedPoint, 8, red, -1);
		}
		*/
	}

	JNIEXPORT void JNICALL Java_com_tum_historicarguide_MainActivity_nativeOpticalDetection(JNIEnv *env, jobject instance, jlong addrRgba, jfloat fps) {
		Mat& mRgb = *(Mat*)addrRgba;

		Mat mGray;
		cvtColor(mRgb, mGray, CV_RGBA2GRAY);

		Mat equalizedMat;
		equalizeHist(mGray, equalizedMat);

		arImgFrames = (int)fps;

		opticalDetection(mRgb, equalizedMat);

		/*
		// Draw touch
		if (touchedPoint.x != 0 || touchedPoint.y != 0) {
			circle(mRgb, touchedPoint, 8, red, -1);
		}
		*/
	}

	// TODO: Not needed anymore
	/*
	JNIEXPORT void JNICALL Java_com_tum_historicarguide_MainActivity_nativeSetTouchPos(JNIEnv *env, jobject instance, jint xCoord, jint yCoord) {
		touchedPoint = Point((int)xCoord, (int)yCoord);
	}
	*/
} // END extern "C" (maybe combining both externs to one)

void generateSearchContours() {
	
	float size = 100;
	vector<Point> tempContour;

	// Generate circle
	float degreePerPoint = 15.0f;
	int numOfPoints = (int)(360.0f / degreePerPoint);
	for (int i = 0; i < numOfPoints; i++) {
		float cosVal = cosf((degreePerPoint * i) * M_PI / 180.0);
		float invCosVal = cosf((90.0f - (degreePerPoint *i)) * M_PI / 180.0);
		tempContour.push_back(Point(center.x - (((size / 2)*cosVal)), center.y - (size / 2)*invCosVal));
	}
	searchContours.push_back(tempContour);

	// Generate Square
	tempContour.clear();
	tempContour.push_back(Point(center.x - (size / 2), center.y - (size / 2)));
	tempContour.push_back(Point(center.x + (size / 2), center.y - (size / 2)));
	tempContour.push_back(Point(center.x + (size / 2), center.y + (size / 2)));
	tempContour.push_back(Point(center.x - (size / 2), center.y + (size / 2)));
	searchContours.push_back(tempContour);

	// Generate rect and parallelograms
	float shearingFactor = 0.25f;
	int numOfParallelograms = (int)(1.0f / shearingFactor);
	for (int i = 0; i < numOfParallelograms; i++) {
		tempContour.clear();
		tempContour.push_back(Point(0, 0));
		tempContour.push_back(Point(size * (shearingFactor*i), size));
		tempContour.push_back(Point(size * (1.5f + (shearingFactor*i)), size));
		tempContour.push_back(Point(size * 1.5f, 0));
		searchContours.push_back(tempContour);
	}

	// Generate arch with extension
	tempContour.clear();
	tempContour = createArch((size / 2), 7, 3, 5);
	searchContours.push_back(tempContour);

	// Generate arch without extension
	tempContour.clear();
	tempContour = createArch((size / 2), 7, 3.0f, 5.0f, 0.33f, false);
	searchContours.push_back(tempContour);

	// Generate semi-circles
	degreePerPoint = 10.0f;
	numOfPoints = ((int)(360.0f / degreePerPoint) / 2) + 1;
	int minPointsInSemiCircle = 5;
	for (int i = minPointsInSemiCircle; i < numOfPoints; i++) {
		tempContour.clear();
		for (int j = 0; j < i; j++) {
			float cosVal = cosf((degreePerPoint * j) * M_PI / 180.0);
			float invCosVal = cosf((90.0f - (degreePerPoint *j)) * M_PI / 180.0);
			tempContour.push_back(Point(center.x - (((size / 2)*cosVal)), center.y - (size / 2)*invCosVal));
		}
		searchContours.push_back(tempContour);
	}
}

void updateSettings(int innerWidth, int innerHeight, int frameWidth, int frameHeight) {

	if (innerFrameWidth != innerWidth) {
		std::ostringstream rowOstr;
		rowOstr << innerWidth;
		rowOstr << " - new height: ";
		rowOstr << innerHeight;
		__android_log_write(ANDROID_LOG_INFO, "new width", rowOstr.str().c_str());
		center = Point((int)innerWidth / 2, (int)innerHeight / 2);
		// touchedPoint = center;

		innerFrameWidth = innerWidth;
		innerFrameHeight = innerHeight;
		cameraFrameWidth = frameWidth;
		cameraFrameHeight = frameHeight;
	}
}

void resetOldMiddleArchProps() {
	old_middle_arch_width = 0;
	old_middle_arch = vector<Point>();
	old_boundingRec = Rect(Point(0, 0), Point(0, 0));
	remainingContourFrames = 0;
}

// Blurs grayscale image and applies sobel filter on it, this results in a binary image (only 1 or 0)
void applySobel(Mat& mGray, int kernel_size) {
	Mat grad_x, grad_y;
	Mat abs_grad_x, abs_grad_y;

	Mat grad;

	int scale = 1;
	int delta = 0;
	int ddepth = CV_16S;

	Mat blurredImg;

	// blur(mGray, blurredImg, Size(kernel_size, kernel_size));
	GaussianBlur(mGray, blurredImg, cv::Size(kernel_size, kernel_size), 3);

	/// Gradient X
	// Scharr(blurredImg, grad_x, ddepth, 0, 1, scale, delta, BORDER_DEFAULT);
	Sobel(blurredImg, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT);
	/// Gradient Y
	// Scharr(blurredImg, grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT);
	// ksize == -1 ==> using Scharr 3x3 Kernel
	Sobel(blurredImg, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT);

	convertScaleAbs(grad_x, abs_grad_x);
	convertScaleAbs(grad_y, abs_grad_y);

	addWeighted(abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad);
	
	Mat kernel = (Mat_<char>(5, 5) << -1, -1, -1, -1, -1,
		-1, 2, 2, 2, -1,
		-1, 2, 8, 2, -1,
		-1, 2, 2, 2, -1,
		-1, -1, -1, -1, -1)/1.55;

	Mat dst;

	filter2D(grad, dst, CV_8U, kernel);

	/*
	Mat image;

	GaussianBlur(grad, image, cv::Size(0, 0), 3);
	addWeighted(grad, 1.25, image, -0.25, 0, image);
	*/
	
	mGray = dst;
}

// Blurs grayscale image and applies canny filter on it, this results in a binary image (only 1 or 0)
void applyBlurAndCanny(Mat& mGray, double factor) {
	Mat detected_edges;

	// Init Config: Size (1,1) + 80, 90

	blur(mGray, detected_edges, Size(3, 3));

	Mat nullMat;

	double otsu_thresh_val = threshold(mGray, nullMat, 0, 255, CV_THRESH_OTSU);

	double high_thresh_val = otsu_thresh_val;
	double lower_thresh_val = otsu_thresh_val * factor;

	Canny(detected_edges, detected_edges, lower_thresh_val, high_thresh_val);	// 60/70


	mGray = detected_edges;
}

void opticalDetectionDebug(Mat& mRgb, Mat& mGray) {
	// mRgb.copyTo(originalRGBImg);

	Mat mGray2;

	mGray.copyTo(mGray2);

	applyBlurAndCanny(mGray, 0.5);

	applySobel(mGray2, 5);

	vector<vector<Point>> firstContours = extractContours(mGray);
	vector<vector<Point>> secondContours = extractContours(mGray2);
	vector<vector<Point>> mergedContours = matchContours(firstContours, secondContours);

	double fixedMatchAccuracy = 0.0375;

	vector<vector<Point>> prunedContours = pruneContours(mRgb, mergedContours, 4, 0.66f, 4.0, true, fixedMatchAccuracy);

	if (!bufferSwitch) {
		switchBuffer[0] = prunedContours;
	}
	else {
		switchBuffer[1] = prunedContours;
	}

	float missingRatio = 0.5f;
	int maxNumOfMisses = (int)floor((arImgFrames * missingRatio) + 0.5f);

	if (!globalContoursMap.empty()) {
		if (!bufferSwitch) {
			globalContoursMap = matchContoursMaps(switchBuffer[0], globalContoursMap, 0.15f, 0.05f, maxNumOfMisses);
		}
		else {
			globalContoursMap = matchContoursMaps(switchBuffer[1], globalContoursMap, 0.15f, 0.05f, maxNumOfMisses);
		}
		vector<vector<Point>> tempContours = contoursMapToContours(globalContoursMap);
		if (tempContours.size() > 0) {
			prevFoundArchContours = tempContours;
			mapARImageOnMiddleArch(mRgb, prevFoundArchContours, archWidth);
			drawFoundContours(mRgb, prevFoundArchContours);
		}
		else if (remainingContourFrames > 0) {
			mapARImageOnMiddleArch(mRgb, prevFoundArchContours, archWidth);
			drawFoundContours(mRgb, prevFoundArchContours);
		} else if (remainingContourFrames <= 0) {
			resetOldMiddleArchProps();
			resizedARImage.release();
		}
	}
	else {
		vector<vector<Point>> combinedContours = matchContours(switchBuffer[0], switchBuffer[1]);
		globalContoursMap = contoursToContoursMap(combinedContours);
		mapARImageOnMiddleArch(mRgb, prevFoundArchContours, archWidth);
		drawFoundContours(mRgb, combinedContours);
	}

	bufferSwitch = !bufferSwitch;

	// drawFoundContours(mRgb, prunedContours);

	drawBorder(mRgb);
}

void opticalDetection(Mat& mRgb, Mat& mGray) {
	mRgb.copyTo(originalRGBImg);

	Mat mGray2;

	mGray.copyTo(mGray2);

	applyBlurAndCanny(mGray, 0.5);

	applySobel(mGray2, 5);

	vector<vector<Point>> firstContours = extractContours(mGray);
	vector<vector<Point>> secondContours = extractContours(mGray2);
	vector<vector<Point>> mergedContours = matchContours(firstContours, secondContours);

	double fixedMatchAccuracy = 0.0375;

	vector<vector<Point>> prunedContours = pruneContours(mRgb, mergedContours, 4, 0.66f, 4.0, true, fixedMatchAccuracy);

	if (!bufferSwitch) {
		switchBuffer[0] = prunedContours;
	}
	else {
		switchBuffer[1] = prunedContours;
	}

	float missingRatio = 0.5f;
	int maxNumOfMisses = (int)floor((arImgFrames * missingRatio) + 0.5f);

	if (!globalContoursMap.empty()) {
		if (!bufferSwitch) {
			globalContoursMap = matchContoursMaps(switchBuffer[0], globalContoursMap, 0.15f, 0.05f, maxNumOfMisses);
		}
		else {
			globalContoursMap = matchContoursMaps(switchBuffer[1], globalContoursMap, 0.15f, 0.05f, maxNumOfMisses);
		}
		vector<vector<Point>> tempContours = contoursMapToContours(globalContoursMap);
		if (tempContours.size() > 0) {
			prevFoundArchContours = tempContours;
			mapARImageOnMiddleArch(mRgb, prevFoundArchContours, archWidth);
		}
		else if (remainingContourFrames > 0) {
			mapARImageOnMiddleArch(mRgb, prevFoundArchContours, archWidth);
		}
		else if (remainingContourFrames <= 0) {
			resetOldMiddleArchProps();
			resizedARImage.release();
		}
	}
	else {
		vector<vector<Point>> combinedContours = matchContours(switchBuffer[0], switchBuffer[1]);
		globalContoursMap = contoursToContoursMap(combinedContours);
		mapARImageOnMiddleArch(mRgb, prevFoundArchContours, archWidth);
	}
	bufferSwitch = !bufferSwitch;
}

Point getCenterOfContour(vector<Point> contour) {
	int sumx = 0;
	int sumy = 0;
	for (int i = 0; i < contour.size(); i++) {
		sumx += contour[i].x;
		sumy += contour[i].y;
	}

	return Point((sumx / contour.size()), sumy / contour.size());
}

float getLongestContourSide(vector<Point> contour) {
	float longestSide = 0;
	for (int i = 0; i < contour.size(); i++) {
		float lengthBetweenPoints = getSideLength(contour[i], contour[(i + 1) % contour.size()]);
		if (lengthBetweenPoints > longestSide) {
			// Check if the longestLine is oriented horizontally
			// Line is horizontally if x values of points differ more than y values
			int xDiff = abs(contour[i].x - contour[(i + 1) % contour.size()].x);
			int yDiff = abs(contour[i].y - contour[(i + 1) % contour.size()].y);
			if (xDiff > yDiff) {
				longestSide = lengthBetweenPoints;
			}
		}
	}
	return longestSide;
}

void mapARImageOnMiddleArch(Mat mRgb, vector<vector<Point>> contours, float middleArchDiameter) {
	vector<Point> middleArch = getAtPointedArch(mRgb, contours);

	// If middle arch was found, map ARImage on Screen
	if (!middleArch.empty()) {
		polylines(mRgb, middleArch, true, blue, 4);
		mapARImageOnZenithOfMiddleArch(mRgb, middleArch, middleArchDiameter, true);
	}
	else if (!old_middle_arch.empty() && remainingContourFrames > 0){
		middleArch = old_middle_arch;
		int tempRemainingContourFrames = remainingContourFrames;
		mapARImageOnZenithOfMiddleArch(mRgb, middleArch, middleArchDiameter, false);
		if (tempRemainingContourFrames <= remainingContourFrames) {
			remainingContourFrames = tempRemainingContourFrames - 1;
		}
	}
	else if (!resizedARImage.empty() && remainingContourFrames > 0){
		resizedARImage.copyTo(glob_mRgb(old_boundingRec));
		remainingContourFrames--;
	}
	else {
		remainingContourFrames--;
	}
	
	if (remainingContourFrames <= 0) {
		resetOldMiddleArchProps();
	}
}

void mapARImageOnZenithOfMiddleArch(Mat mRgb, vector<Point> middleArch, float middleArchDiameter, bool newArchFound) {
	// Transform Point to destination Mat
	if (newArchFound) {
		middleArch = convertContourToOtherMat(mRgb, glob_mRgb, middleArch);
	}
	Point zenithOfArch = getZenithOfContour(middleArch);
	int diameterOfFoundArch = boundingRect(middleArch).width;
	// Scale ARImage
	int rectWidth = (int)(diameterOfFoundArch * ((float)arImage.cols / middleArchDiameter));
	float imgRatio = (float)arImage.rows / (float)arImage.cols;
	int rectHeight = (int)((float)rectWidth * imgRatio);
	Size size(rectWidth, rectHeight);
	Mat preScaledImg;
	resize(arImage, preScaledImg, size);
	// Get Rect to place ARImage in
	Point topLeft = Point(max(zenithOfArch.x - (int)(rectWidth * archZenithPercPosX), 0), max(zenithOfArch.y - (int)(rectHeight * archZenithPercPosY), 0));
	Point bottomRight = Point(min(zenithOfArch.x + (int)(rectWidth * (1.0f - archZenithPercPosX)), cameraFrameWidth), (min(zenithOfArch.y + (int)(rectHeight * (1.0f - archZenithPercPosY)), cameraFrameHeight)));
	Rect boundingRec = Rect(topLeft, bottomRight);
	// Crop ARImage if it would be out of bounds, because of translation or big scaling
	Point topLeftDiffAR = Point(min(zenithOfArch.x - (int)(rectWidth * archZenithPercPosX), 0), min(zenithOfArch.y - (int)(rectHeight * archZenithPercPosY), 0));
	Point topLeftAR = Point(-topLeftDiffAR.x, -topLeftDiffAR.y);
	Rect preCroppedRec = Rect(topLeftAR, Size(boundingRec.width, boundingRec.height));
	Mat preCroppedImg = preScaledImg(preCroppedRec);
	// Change Alpha Value of ARImage and copy it to RGB frame
	Mat backgroundImg = originalRGBImg(boundingRec);
	setBlackPixelsTransparent(preCroppedImg, backgroundImg);
	changeAlphaValue(preCroppedImg, backgroundImg, 0.9f);
	resizedARImage = preCroppedImg;
	old_boundingRec = boundingRec;
	resizedARImage.copyTo(glob_mRgb(boundingRec));
	old_middle_arch = middleArch;
	remainingContourFrames = arImgFrames * minSecsToShowARImg;
	polylines(glob_mRgb, middleArch, true, red, 8);
}

vector<Point> convertContourToOtherMat(Mat srcMat, Mat dstMat, vector<Point> contour) {
	Point dstMatCenter = Point(dstMat.cols/2, dstMat.rows/2);
	Point topLeftOfInnerFrame = Point(dstMatCenter.x - (innerFrameWidth/2), dstMatCenter.y - (innerFrameHeight/2));
	circle(dstMat, dstMatCenter, 8, blue, -1);
	circle(dstMat, topLeftOfInnerFrame, 8, green, -1);
	vector<Point> convertedContour;
	for (Point p : contour) {
		Point convertedPoint = Point(p.x + topLeftOfInnerFrame.x, p.y + topLeftOfInnerFrame.y);
		convertedContour.push_back(convertedPoint);
	}
	return convertedContour;
}

Point getZenithOfContour(vector<Point> contour) {
	Point highestPoint = Point(0, cameraFrameHeight);
	for (Point p : contour) {
		if (p.y < highestPoint.y) {
			highestPoint = p;
		}
	}
	Point contourCenter = getCenterOfContour(contour);
	return Point(contourCenter.x, highestPoint.y);
}

vector<Point> getAtPointedArch(Mat mRgb, vector<vector<Point>> contours) {
	// Check if onScreen was tapped

	// TODO: Wird eigentlich im neuen Ansatz nicht mehr benötigt, da alle Konturen im inneren frame in frage kommen
	/*
	vector<vector<Point>> contoursPointedAt;
	if (touchedPoint.x != 0 || touchedPoint.y != 0) {
		// Get all contours that were hit by the touchedPoint
		for (vector<Point> contour : contours) {
			Rect boundingRec = boundingRect(contour);
			// extend boundingRec to increase detecting area
			Point topLeft = Point(max(boundingRec.tl().x - 5, 0), max(boundingRec.tl().y - 75, 0));
			Point bottomRight = Point(min(boundingRec.br().x + 5, innerFrameWidth), min(boundingRec.br().y + 75, innerFrameHeight));
			boundingRec = Rect(topLeft, bottomRight);
			if (pointInRect(touchedPoint, boundingRec)) {
				contoursPointedAt.push_back(contour);
			}
		}
	} else {
		return vector<Point>();
	}
	*/
	vector<Point> biggestContour;
	int maxWidth = -1;
	int widthSum = 0;
	// Get biggest contour
	for (vector<Point> contour : contours) {
		Rect boundingRec = boundingRect(contour);
		int currentWidth = boundingRec.width;
		if (contourIsArch(contour)) {
			widthSum += currentWidth;
			// Current contour is bigger than the biggest found contour
			if (currentWidth > maxWidth) {
				// Contour is at least 20% of innerFrameWidth
				if (boundingRec.width > innerFrameWidth * 0.2) {
					if (boundingRec.width * 0.3 < boundingRec.height) {
						biggestContour = contour;
						maxWidth = currentWidth;
					}
				}
			}
		}
	}

	vector<vector<Point>> potentialContours;
	vector<Point> mostCenterContour;
	double minDistToCenter = innerFrameWidth;
	// Get contours with size similar to biggest contour
	if (!biggestContour.empty()) {
		int minWidth = maxWidth * 0.8;
		int potentialContoursWidthSum = 0;
		for (vector<Point> contour : contours) {
			Rect boundingRec = boundingRect(contour);
			if (boundingRec.width > minWidth && contourIsArch(contour)) {
				potentialContours.push_back(contour);
				potentialContoursWidthSum += boundingRec.width;
				// Get potential contour that is closest to middle
				Point centerOfContour = getCenterOfContour(contour);
				double distanceToCenter = norm(centerOfContour - center);
				if (distanceToCenter < minDistToCenter) {
					minDistToCenter = distanceToCenter;
					mostCenterContour = contour;
				}
			}
		}
		float averageWidth = (float)potentialContoursWidthSum / (float)potentialContours.size();
		last_arch_widths[last_pushed_index] = averageWidth;
	}
	else {
		last_arch_widths[last_pushed_index] = 0;
	}
	last_pushed_index = (last_pushed_index + 1) % number_of_arch_widths;


	// Hier nehmen wir die Summe der gefundenen Breiten
	widthSum = 0;
	int numberOfWidths = 0;
	for (int width : last_arch_widths) {
		if (width != 0) {
			widthSum += width;
			numberOfWidths++;
		}
	}
	if (numberOfWidths > 2) {
		float avWidth = (float)widthSum / (float)numberOfWidths;
		// Get the highets/lowest contour depending if the arch has a background or not
		vector<Point> foundContour;
		if (archHasBackground) {
			// e.g. Feldherrnhalle ==> return contour with highest zenith
			Point highestZenith = Point(innerFrameWidth, innerFrameHeight);
			for (vector<Point> potentialContour : potentialContours) {
				Point currentZenith = getZenithOfContour(potentialContour);
				if (currentZenith.y < highestZenith.y) {
					Point centerOfContour = getCenterOfContour(potentialContour);
					double distanceToCenter = norm(centerOfContour - center);
					// Check if the found contour is close to the center
					if (distanceToCenter <= minDistToCenter * 1.25) {
						highestZenith = currentZenith;
						foundContour = potentialContour;
					}
				}
			}
		}
		else {
			// e.g. Siegestor ==> return contour with lowest zenith
			Point lowestZenith = Point(0, 0);
			for (vector<Point> potentialContour : potentialContours) {
				Point currentZenith = getZenithOfContour(potentialContour);
				if (currentZenith.y > lowestZenith.y) {
					Point centerOfContour = getCenterOfContour(potentialContour);
					double distanceToCenter = norm(centerOfContour - center);
					// Check if the found contour is close to the center
					if (distanceToCenter <= minDistToCenter * 1.25) {
						lowestZenith = currentZenith;
						foundContour = potentialContour;
					}
				}
			}
		}

		// It's ensured, that foundContour's width is at least, 80 percent of the biggestContour
		// Its's also ensured, that the highest contour will be found if the arch has a background and the lowest if it has no background
		if (!foundContour.empty()) {
			int foundContWidth = boundingRect(foundContour).width;
			// Checks if found contour is similar to previous found arch widths and if found contour isn't very small
			if (avWidth * 0.9 < foundContWidth && avWidth * 1.1 > foundContWidth && innerFrameWidth * 0.25 < foundContWidth) {
				return biggestContour;
			}
			else {
				return vector<Point>();
			}
		}

		// Return biggestContour if it has been found
		if (!biggestContour.empty()) {
			int biggestContWidth = boundingRect(biggestContour).width;
			// Checks if biggest contour is similar to previous found arch widths and if biggest contour isn't very small
			if (avWidth * 0.9 < biggestContWidth && avWidth * 1.1 > biggestContWidth && innerFrameWidth * 0.25 < biggestContWidth) {
				return biggestContour;
			}
			else {
				return vector<Point>();
			}
		}
		else {
			return vector<Point>();
		}
	}
	else {
		return vector<Point>();
	}
}

void setBlackPixelsTransparent(Mat image, Mat backgroundImg) {
	for (int i = 0; i < image.cols; i++) {
		for (int j = 0; j < image.rows; j++) {
			// Search for green pixels
			if (image.at<cv::Vec4b>(j, i)[0] <= 5 && image.at<cv::Vec4b>(j, i)[1] <= 5 && image.at<cv::Vec4b>(j, i)[2] <= 5) {
				image.at<cv::Vec4b>(j, i)[0] = backgroundImg.at<cv::Vec4b>(j, i)[0];
				image.at<cv::Vec4b>(j, i)[1] = backgroundImg.at<cv::Vec4b>(j, i)[1];
				image.at<cv::Vec4b>(j, i)[2] = backgroundImg.at<cv::Vec4b>(j, i)[2];
			}
		}
	}
}

void changeAlphaValue(Mat image, Mat backgroundImg, float opaquePercentage) {
	// only scale x percent
	float eachSidePercentage = (1.0f - opaquePercentage) / 2.0f;
	int borderStrength = (int)(min(image.rows, image.cols) * eachSidePercentage);
	Point topLeft = Point(borderStrength, borderStrength);
	Point topRight = Point(image.cols - borderStrength, borderStrength);
	Point bottomLeft = Point(borderStrength, image.rows - borderStrength);
	Point bottomRight = Point(image.cols - borderStrength, image.rows - borderStrength);
	if (image.channels() == 4 && opaquePercentage < 1.0f && opaquePercentage > 0.0f) {
		// 3 cases needed, only top and bottom, only left and right, only corners
		for (int i = 0; i < image.cols; i++) {
			if (i < borderStrength) {
				for (int j = 0; j < image.rows; j++) {
					// only corner top left
					if (j < borderStrength) {
						float tempC = sqrt(pow((float)i - topLeft.x, 2) + pow((float)j - topLeft.y, 2));
						float perc;
						if (tempC >= borderStrength) {
							perc = 1;
						} else {
							perc = tempC / borderStrength;
						}
						image.at<cv::Vec4b>(j, i)[0] = perc * backgroundImg.at<cv::Vec4b>(j, i)[0] + (1 - perc) * image.at<cv::Vec4b>(j, i)[0];
						image.at<cv::Vec4b>(j, i)[1] = perc * backgroundImg.at<cv::Vec4b>(j, i)[1] + (1 - perc) * image.at<cv::Vec4b>(j, i)[1];
						image.at<cv::Vec4b>(j, i)[2] = perc * backgroundImg.at<cv::Vec4b>(j, i)[2] + (1 - perc) * image.at<cv::Vec4b>(j, i)[2];
					} 
					// only corner bottom left
					else if (j > image.rows - borderStrength) {
						float tempC = sqrt(pow((float)i - bottomLeft.x, 2) + pow((float)j - bottomLeft.y, 2));
						float perc;
						if (tempC >= borderStrength) {
							perc = 1;
						}
						else {
							perc = tempC / borderStrength;
						}
						image.at<cv::Vec4b>(j, i)[0] = perc * backgroundImg.at<cv::Vec4b>(j, i)[0] + (1 - perc) * image.at<cv::Vec4b>(j, i)[0];
						image.at<cv::Vec4b>(j, i)[1] = perc * backgroundImg.at<cv::Vec4b>(j, i)[1] + (1 - perc) * image.at<cv::Vec4b>(j, i)[1];
						image.at<cv::Vec4b>(j, i)[2] = perc * backgroundImg.at<cv::Vec4b>(j, i)[2] + (1 - perc) * image.at<cv::Vec4b>(j, i)[2];
					}
					// only left side
					else {
						float perc = (float)i / (float)borderStrength;
						image.at<cv::Vec4b>(j, i)[0] = (1 - perc) * backgroundImg.at<cv::Vec4b>(j, i)[0] + perc * image.at<cv::Vec4b>(j, i)[0];
						image.at<cv::Vec4b>(j, i)[1] = (1 - perc) * backgroundImg.at<cv::Vec4b>(j, i)[1] + perc * image.at<cv::Vec4b>(j, i)[1];
						image.at<cv::Vec4b>(j, i)[2] = (1 - perc) * backgroundImg.at<cv::Vec4b>(j, i)[2] + perc * image.at<cv::Vec4b>(j, i)[2];
					}
				}
			}
			else if (i > image.cols - borderStrength) {
				for (int j = 0; j < image.rows; j++) {
					// only corner top right
					if (j < borderStrength) {
						float tempC = sqrt(pow((float)i - topRight.x, 2) + pow((float)j - topRight.y, 2));
						float perc;
						if (tempC >= borderStrength) {
							perc = 1;
						}
						else {
							perc = tempC / borderStrength;
						}
						image.at<cv::Vec4b>(j, i)[0] = perc * backgroundImg.at<cv::Vec4b>(j, i)[0] + (1 - perc) * image.at<cv::Vec4b>(j, i)[0];
						image.at<cv::Vec4b>(j, i)[1] = perc * backgroundImg.at<cv::Vec4b>(j, i)[1] + (1 - perc) * image.at<cv::Vec4b>(j, i)[1];
						image.at<cv::Vec4b>(j, i)[2] = perc * backgroundImg.at<cv::Vec4b>(j, i)[2] + (1 - perc) * image.at<cv::Vec4b>(j, i)[2];
					} 
					// only corner bottom right 
					else if (j > image.rows - borderStrength) {
						float tempC = sqrt(pow((float)i - bottomRight.x, 2) + pow((float)j - bottomRight.y, 2));
						float perc;
						if (tempC >= borderStrength) {
							perc = 1;
						}
						else {
							perc = tempC / borderStrength;
						}
						image.at<cv::Vec4b>(j, i)[0] = perc * backgroundImg.at<cv::Vec4b>(j, i)[0] + (1 - perc) * image.at<cv::Vec4b>(j, i)[0];
						image.at<cv::Vec4b>(j, i)[1] = perc * backgroundImg.at<cv::Vec4b>(j, i)[1] + (1 - perc) * image.at<cv::Vec4b>(j, i)[1];
						image.at<cv::Vec4b>(j, i)[2] = perc * backgroundImg.at<cv::Vec4b>(j, i)[2] + (1 - perc) * image.at<cv::Vec4b>(j, i)[2];
					}
					// only right side
					else {
						float perc = (float)(image.cols - i) / (float)borderStrength;
						image.at<cv::Vec4b>(j, i)[0] = (1 - perc) * backgroundImg.at<cv::Vec4b>(j, i)[0] + perc * image.at<cv::Vec4b>(j, i)[0];
						image.at<cv::Vec4b>(j, i)[1] = (1 - perc) * backgroundImg.at<cv::Vec4b>(j, i)[1] + perc * image.at<cv::Vec4b>(j, i)[1];
						image.at<cv::Vec4b>(j, i)[2] = (1 - perc) * backgroundImg.at<cv::Vec4b>(j, i)[2] + perc * image.at<cv::Vec4b>(j, i)[2];
					}
				}
			}
			else {
				for (int j = 0; j < image.rows; j++) {
					// only top
					if (j < borderStrength) {
						float perc = (float)j / (float)borderStrength;
						image.at<cv::Vec4b>(j, i)[0] = (1 - perc) * backgroundImg.at<cv::Vec4b>(j, i)[0] + perc * image.at<cv::Vec4b>(j, i)[0];
						image.at<cv::Vec4b>(j, i)[1] = (1 - perc) * backgroundImg.at<cv::Vec4b>(j, i)[1] + perc * image.at<cv::Vec4b>(j, i)[1];
						image.at<cv::Vec4b>(j, i)[2] = (1 - perc) * backgroundImg.at<cv::Vec4b>(j, i)[2] + perc * image.at<cv::Vec4b>(j, i)[2];
					} 
					// only bottom
					else if (j > image.rows - borderStrength) {
						float perc = (float)(image.rows - j) / (float)borderStrength;
						image.at<cv::Vec4b>(j, i)[0] = (1 - perc) * backgroundImg.at<cv::Vec4b>(j, i)[0] + perc * image.at<cv::Vec4b>(j, i)[0];
						image.at<cv::Vec4b>(j, i)[1] = (1 - perc) * backgroundImg.at<cv::Vec4b>(j, i)[1] + perc * image.at<cv::Vec4b>(j, i)[1];
						image.at<cv::Vec4b>(j, i)[2] = (1 - perc) * backgroundImg.at<cv::Vec4b>(j, i)[2] + perc * image.at<cv::Vec4b>(j, i)[2];
					}
				}
			}
		}
	}
}

// Create border around frame to be analyzed
void drawBorder(Mat mRgb, int borderSize) {
	Point leftTop = Point(borderSize / 2, borderSize / 2);
	Point rightBottom = Point(mRgb.cols - borderSize / 2, mRgb.rows - borderSize / 2);
	rectangle(mRgb, leftTop, rightBottom, white, borderSize);
}

vector<Point> createArch(int radius, int numOfPointsInSemiCircle, float scaleX, float scaleY, float circleScaleX, bool archExtension) {
	vector<Point> arch;

	
	if (archExtension) {
		// Right extension
		arch.push_back(Point(center.x + ((radius / 2)*scaleX), center.y + ((radius / 2)*scaleY)));
		// Left extension
		arch.push_back(Point(center.x - ((radius / 2)*scaleX), center.y + ((radius / 2)*scaleY)));
	}	

	if (numOfPointsInSemiCircle > 4) {
		float degreePerPoint = 180.0f / (numOfPointsInSemiCircle - 1);
		// std::ostringstream rowOstr2;
		// rowOstr2 << degreePerPoint;
		// __android_log_write(ANDROID_LOG_INFO, "degreePerPoint", rowOstr2.str().c_str());
		for (int i = 0; i < numOfPointsInSemiCircle; i++) {
			float cosVal = cosf((degreePerPoint * i) * M_PI / 180.0);
			float invCosVal = cosf((90.0f - (degreePerPoint *i)) * M_PI / 180.0);
			if (i == 0 || i == numOfPointsInSemiCircle - 1) {
				arch.push_back(Point(center.x - (((radius / 2)*cosVal)*scaleX), center.y - (((radius / 2)*invCosVal)*scaleY)));
			} else {
				arch.push_back(Point(center.x - (((radius / 2)*cosVal)*(scaleX + circleScaleX)), center.y - (((radius / 2)*invCosVal)*scaleY)));
			}
		}
	}
	else {
		// Create Arch with 5 points
		float degreePerPoint = 180.0f / 4;
		for (int i = 0; i < 5; i++) {
			arch.push_back(Point(center.x - (((radius / 2)*cosf((degreePerPoint * i) * M_PI / 180.0))*scaleX), center.y - (((radius / 2)*cosf((degreePerPoint * (numOfPointsInSemiCircle - i)) * M_PI / 180.0))*scaleY)));
		}
	}
	return arch;
}

vector<vector<Point>> eliminateContoursOutOfBounds(vector<vector<Point>> contours, int width, int heigth) {
	vector<vector<Point>> adjustedContours;
	int maxNumOfMisses = contours.size() / 2;
	for (vector<Point> contour : contours) {
		int counter = 0;
		for (Point point : contour) {
			if (point.x < 0 || point.x > width) {
				counter++;
			} else if (point.y < 0 || point.y > heigth) {
				counter++;
			}
		}
		if (counter < maxNumOfMisses) {
			adjustedContours.push_back(contour);
		}
	}
	return adjustedContours;
}

void drawFoundContours(Mat srcDst, vector<vector<Point>> foundContours, int thickness) {
	for (vector<Point> foundContour : foundContours) {
		Scalar shapeColor;
		int numOfPointsInContour = foundContour.size();
		switch (numOfPointsInContour)
		{
		case 3:
			shapeColor = Scalar(0, 0, 255);
			break;
		case 4:
			shapeColor = Scalar(255, 0, 0);
			break;
		case 5:
			shapeColor = Scalar(255, 0, 255);
			break;
		default:
			shapeColor = Scalar(255, 255, 255);
			break;
		}
		if (contourIsArch(foundContour)) {
			for (int i = 5; i < searchContours.size(); i++) {
				if (matchShapes(foundContour, searchContours[i], CV_CONTOURS_MATCH_I1, 0) < 0.15) {
					shapeColor = Scalar(77, 255, 0);
				}
			}
		} else if (matchShapes(foundContour, searchContours[0], CV_CONTOURS_MATCH_I1, 0) < 0.03) {
			shapeColor = Scalar(255, 77, 127);
		}
		Point zenith = getZenithOfContour(foundContour);
		circle(srcDst, zenith, 8, red, -1);
		polylines(srcDst, foundContour, true, shapeColor, thickness);
		// drawContours(mRgb, mergedContours, i, shapeColor, 2, 8, vector<Vec4i>(), 0, Point());
	}
}

vector<vector<Point>> extractArches(Mat mRgb, vector<vector<Point>> contours) {
	vector<vector<Point>> extractedContours;
	for (vector<Point> contour : contours) {
		if (contour.size() > 4) {
			if (contourIsArch(contour)) {
				extractedContours.push_back(contour);
			}
		}
	}
	return extractedContours;
}

bool contourIsArch(vector<Point> contour, bool drawPoints) {
	// Find sharp angles
	int numOfPoints = contour.size();
	if (numOfPoints < 6) {
		return false;
	}
	vector<int> acuteIndices;
	int leftAcuteAngles = 2;
	for (int i = 0; i < numOfPoints; i++) {
		// acosf(0.9732f)  * 180.0 / M_PI
		float sideA = getSideLength(contour[(i - 1) % numOfPoints], contour[i]);
		float sideB = getSideLength(contour[i], contour[(i + 1) % numOfPoints]);
		float sideC = getSideLength(contour[(i - 1) % numOfPoints], contour[(i + 1) % numOfPoints]);
		float angle = getOppositeAngleOfSideC(sideA, sideB, sideC);
		/*
		if (angle < 120 && angle > 60) {
		if (leftAcuteAngles > 0) {
		leftAcuteAngles--;
		acuteIndices.push_back(i);
		}
		else {
		return false;
		}
		}
		*/

		if (angle > 115) {
			if (drawPoints) {
				// circle(mRgb, contour[i], 8, Scalar(0, 247, 255), 2, 8, 0);
			}
		}
		else {
			// Acute/Sharp Angle
			if (drawPoints) {
				// circle(mRgb, contour[i], 8, Scalar(255, 179, 0), 2, 8, 0);
			}
			if (leftAcuteAngles > 0) {
				leftAcuteAngles--;
				acuteIndices.push_back(i);
			}
			else {
				return false;
			}
			
		}
	}
	if (acuteIndices.size() == 2) {
		// Check if indices are next to each other
		if (abs(acuteIndices[0] - acuteIndices[1]) == 1) {
			if (pointIsAtBottom(contour[acuteIndices[0]], contour) && pointIsAtBottom(contour[acuteIndices[1]], contour)) {
				if (pointsAreAtOuterPos(contour[acuteIndices[0]], contour[acuteIndices[1]], contour)) {
					// Draw circle for each point
					for (int i = 0; i < contour.size(); i++) {
						/*
						if (i == acuteIndices[0] || i == acuteIndices[1]) {
							circle(mRgb, contour[i], 8, Scalar(255, 179, 0), 2, 8, 0);
						}
						else {
							circle(mRgb, contour[i], 8, Scalar(0, 247, 255), 2, 8, 0);
						}
						*/
					}
					return true;
				}
				else {
					return false;
				}
			} else {
				return false;
			}
		}
		else {
			if (abs(((acuteIndices[0] + 1) % numOfPoints) - ((acuteIndices[1] + 1) % numOfPoints)) == 1) {
				if (pointIsAtBottom(contour[((acuteIndices[0] + 1) % numOfPoints)], contour) && pointIsAtBottom(contour[((acuteIndices[1] + 1) % numOfPoints)], contour)) {
					if (pointsAreAtOuterPos(contour[((acuteIndices[0] + 1) % numOfPoints)], contour[((acuteIndices[1] + 1) % numOfPoints)], contour)) {
						// Draw circle for each point
						for (int i = 0; i < contour.size(); i++) {
							/*
							if (i == acuteIndices[0] || i == acuteIndices[1]) {
								circle(mRgb, contour[i], 8, Scalar(255, 179, 0), 2, 8, 0);
							}
							else {
								circle(mRgb, contour[i], 8, Scalar(0, 247, 255), 2, 8, 0);
							}
							*/
						}
						return true;
					}
					else {
						return false;
					}
				}
				else {
					return false;
				}
			}
			else {
				return false;
			}
		}
	}
	else {
		return false;
	}
}

bool pointsAreAtOuterPos(Point point1, Point point2, vector<Point> contour) {
	int counter = 1;
	Point leftPoint;
	Point rightPoint;
	if (point1.x < point2.x) {
		leftPoint = point1;
		rightPoint = point2;
	}
	else {
		leftPoint = point2;
		rightPoint = point1;
	}

	for (int i = 0; i < contour.size(); i++) {
		int variance = 5;
		if (leftPoint.x > contour[i].x + variance) {
			// If a point is closer to the left border than the leftmost-Point decrease counter
			counter--;
		}
		if (rightPoint.x < contour[i].x - variance) {
			// If a point is closer to the right border than the rightmost-Point decrease counter
			counter--;
		}
	}
	if (counter < 0) {
		return false;
	}
	else {
		return true;
	}
}

bool pointIsAtBottom(Point point, vector<Point> contour) {
	int counter = 1;
	for (int i = 0; i < contour.size(); i++) {
		if (point.y < contour[i].y) {
			// If a point is found which is closer to the bottom than the checked point, decrease the counter
			counter--;
		}
	}
	if (counter < 0) {
		return false;
	}
	else {
		return true;
	}
}

float getSideLength(Point p1, Point p2) {
	float xDiff = p1.x - p2.x;
	float yDiff = p1.y - p2.y;
	float sideLength = sqrt(pow(xDiff, 2) + pow(yDiff, 2));
	return sideLength;
}

// Calculates the angle in the triangle using the law of cosines for sideC
// c^2 = a^2 + b^2 - 2*a*b*cos(gamma)
// <==> c^2 - a^2 - b^2 = - 2*a*b*cos(gamma)
// <==> -(c^2 - a^2 - b^2) = 2*a*b*cos(gamma)
// <==> powDiff = 2*a*b*cos(gamma)
// <==> powDiff = multSides*cos(gamma)
// <==> powDiff/multSides = cos(gamma)
// <==> gamma = acos(powDiff/multSides)
// Returns opposite angle of side c
// TODO: Müsste eigentlich "lawOfCosines" heißen
float getOppositeAngleOfSideC(float sideA, float sideB, float sideC) {
	float powDiff = -(pow(sideC, 2) - pow(sideA, 2) - pow(sideB, 2));
	float multSides = 2 * sideA * sideB;
	float angle = acosf(powDiff / multSides) * 180.0 / M_PI;
	return angle;
}

bool contoursRoughlyAtSamePos(vector<Point> contour1, vector<Point> contour2, float maxDeviation, int maxTotalDiffer) {
	int remainingMisses = maxTotalDiffer;
	if (abs(contour1.size() - contour2.size()) < maxTotalDiffer) {
		for (int i = 0; i < contour1.size(); i++) {
			for (int j = 0; j < contour2.size(); j++) {
				// If a point at roughly same pos is found exit the inner for loop
				if (pointsAtSamePos(contour1[i], contour2[j], maxDeviation)) {
					break;
				} 
				// If the point couldn't be found decrease the counter
				else if (j == contour2.size()-1) {
					// std::ostringstream rowOstr11;
					// rowOstr11 << 0;
					// __android_log_write(ANDROID_LOG_INFO, "FAIL", rowOstr11.str().c_str());
					remainingMisses--;
				}
			}
		}
		/*
		std::ostringstream rowOstr4;
		rowOstr4 << remainingMisses;
		__android_log_write(ANDROID_LOG_INFO, "remainingMisses", rowOstr4.str().c_str());
		*/
		if (remainingMisses >= 0) {
			return true;
		}
		else
		{
			return false;
		}
	}
	else 
	{
		return false;
	}
}

// Simple map not possible, because unique keys are needed 
// The remainingNumberOfMisses is not unique so can't be used as map key
// The contour is possibly unique but you need to introduce operators which dont work
ContoursMap contoursToContoursMap(vector<vector<Point>> contours, int maxNumOfMisses) {
	ContoursMap contoursMap;
	for (int i = 0; i < contours.size(); i++) {
		contoursMap[i] = make_pair(contours[i], maxNumOfMisses);
	}
	return contoursMap;
}

// Get contours in ContoursMap
vector<vector<Point>> contoursMapToContours(ContoursMap contoursMap) {
	vector<vector<Point>> contours;
	ContoursMap::iterator it = contoursMap.begin();
	for (it = contoursMap.begin(); it != contoursMap.end(); it++) {
		if (it->second.second > 0) {
			contours.push_back(it->second.first);
		}
	}
	return contours;
}

ContoursMap matchContoursMaps(vector<vector<Point>> newContours, ContoursMap oldContoursMap, float maxDeviation, float matchShapeDeviation, int maxNumOfMisses) {
	ContoursMap combinedContoursMap;
	int val = 0;
	ContoursMap::iterator it = oldContoursMap.begin();
	for (it = oldContoursMap.begin(); it != oldContoursMap.end(); it++) {
		bool alreadyAdded = false;
		for (int j = 0; j < newContours.size(); j++) {
			if (contourCentersMatching(newContours[j], it->second.first)) {
				// Contour j has match with contoursMap contour i
				combinedContoursMap[val] = make_pair(newContours[j], min(it->second.second + 1, maxNumOfMisses));
				newContours.erase(newContours.begin() + j);
				alreadyAdded = true;
				val++;
				break;
			}
		}
		// contoursMap i couldn't be matched
		if (!alreadyAdded) {
			it->second.second--;
			// Add old contour if it's still "alive"
			if (it->second.second > 0) {
				combinedContoursMap[val] = make_pair(it->second.first, it->second.second);
				val++;
			}
		}
	}

	int size = combinedContoursMap.size();
	// Add all new Contours which couldn't be matched
	for (int i = 0; i < newContours.size(); i++) {
		combinedContoursMap[size + i] = make_pair(newContours[i], maxNumOfMisses/2);
	}

	oldContoursMap.clear();

	return combinedContoursMap;
}

bool contourCentersMatching(vector<Point> contour1, vector<Point> contour2, float maxBoundAreaDiff, float maxCenterDiff) {
	Rect contour1Rec = boundingRect(contour1);
	Rect contour2Rec = boundingRect(contour2);
	int smallerArea = min(contour1Rec.area(), contour2Rec.area());
	int biggerArea = max(contour1Rec.area(), contour2Rec.area());
	// Check if boundingArea is of equal Size
	if ((float)(smallerArea / biggerArea) > maxBoundAreaDiff) {
		// Check if the centers are at similar Pos
		Point smallerContourCenter;
		Point biggerContourCenter;
		if (smallerArea == contour1Rec.area()) {
			// contour1 < contour2
			smallerContourCenter = getCenterOfContour(contour1);
			biggerContourCenter = getCenterOfContour(contour2);
			Point topLeft = Point(biggerContourCenter.x - contour2Rec.width * maxCenterDiff, biggerContourCenter.y - contour2Rec.height * maxCenterDiff);
			Point bottomRight = Point(biggerContourCenter.x + contour2Rec.width * maxCenterDiff, biggerContourCenter.y + contour2Rec.height * maxCenterDiff);
			Rect innerRec = Rect(topLeft, bottomRight);
			if (pointInRect(smallerContourCenter, innerRec)) {
				return true;
			} else {
				return false;
			}
		}
		else {
			// contour2 < contour1
			smallerContourCenter = getCenterOfContour(contour2);
			biggerContourCenter = getCenterOfContour(contour1);
			Point topLeft = Point(biggerContourCenter.x - contour1Rec.width * maxCenterDiff, biggerContourCenter.y - contour1Rec.height * maxCenterDiff);
			Point bottomRight = Point(biggerContourCenter.x + contour1Rec.width * maxCenterDiff, biggerContourCenter.y + contour1Rec.height * maxCenterDiff);
			Rect innerRec = Rect(topLeft, bottomRight);
			if (pointInRect(smallerContourCenter, innerRec)) {
				return true;
			}
			else {
				return false;
			}
		}
	}
	else {
		return false;
	}
}

bool pointInRect(Point point, Rect rect) {
	if (point.x >= rect.tl().x && point.x <= rect.br().x) {
		if (point.y >= rect.tl().y && point.y <= rect.br().y) {
			return true;
		}
		else {
			return false;
		}
	} else {
		return false;
	}
}

vector<vector<Point>> eliminateDuplicateContours(vector<vector<Point>> contours) {
	vector<vector<Point>> copiedContours = contours;
	for (int i = 0; i < contours.size() - 1; i++) {
		for (int j = i + 1; j < copiedContours.size(); j++) {
			if (matchShapes(contours[i], copiedContours[j], CV_CONTOURS_MATCH_I1, 0) < 0.08f) {
				if (contoursRoughlyAtSamePos(contours[i], copiedContours[j], 0.25f)) {
					copiedContours.erase(copiedContours.begin() + j);
					break;
				}
			}
		}
	}
	return copiedContours;
}

vector<vector<Point>> matchContours(vector<vector<Point>> contours1, vector<vector<Point>> contours2, float maxDeviation, float matchShapeDeviation) {
	vector<vector<Point>> combinedContours;

	for (int i = 0; i < contours1.size(); i++) {
		combinedContours.push_back(contours1[i]);
	}
	for (int i = 0; i < contours2.size(); i++) {
		combinedContours.push_back(contours2[i]);
	}
	if (!combinedContours.empty()) {
		combinedContours = eliminateDuplicateContours(combinedContours);
	}

	return combinedContours;
}

vector<vector<Point>> extractContours(Mat grayImg, int minShapeArea, int minSideLength, int minVerts, int maxVerts) {
	vector<vector<Point> > contours;
	vector<Vec4i> hierarchy;

	findContours(grayImg, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

	// Create convex hull and approximate contours
	vector<vector<Point> > contours_poly(contours.size());
	vector<vector<Point>> hull(contours.size());
	vector<vector<Point>> extracted_contours(0);

	// Filtering very small contours and drawing the remaining contours (+ pushing them to another vector)
	for (int i = 0; i < contours.size(); i++) {
		if (contourArea(contours[i]) > minShapeArea) {
			convexHull(Mat(contours[i]), hull[i], false);
			approxPolyDP(Mat(hull[i]), contours_poly[i], 2, true);
			int numOfPointsInContour = contours_poly[i].size();
			Rect rect = boundingRect(contours_poly[i]);
			if (rect.height > minSideLength && rect.width > minSideLength && numOfPointsInContour <= maxVerts && numOfPointsInContour >= minVerts) {
				extracted_contours.push_back(contours_poly[i]);
			}
		}
	}
	return extracted_contours;
}

vector<vector<Point>> pruneContours(Mat mRgb, vector<vector<Point>> contours, int minPoints, float maxRotArea, double maxProportion, bool onlySpecContours, double matchAccuracy) {
	vector<vector<Point>> prunedContours;

	for (vector<Point> contour : contours) {
		RotatedRect rotatedRect = minAreaRect(contour);
		Point2f vertices[4];
		rotatedRect.points(vertices);
		double dist1 = norm(vertices[0] - vertices[1]);
		double dist2 = norm(vertices[1] - vertices[2]);

		// Find top/bottom side and left/right side
		Rect boundRec = boundingRect(contour);
		int rotRecArea = rotatedRect.size.area();
		int bounRecArea = boundRec.area();
		if (contour.size() >= minPoints) {
			if (((float)rotRecArea / (float)bounRecArea) > maxRotArea) {
				if (max(dist1, dist2) / min(dist1, dist2) < maxProportion) {
					if (onlySpecContours) {
						bool isSimilarShape = false;
						// Test if contour has similar shape than one of the first three searchContours
						for (int i = 0; i < 3; i++) {
							Point center = getCenterOfContour(contour);
							if (matchShapes(contour, searchContours[i], CV_CONTOURS_MATCH_I1, 0) < matchAccuracy) {
								isSimilarShape = true;
								break;
							}
						}
						// Add contour if it has desired shape
						if (isSimilarShape) {
							prunedContours.push_back(contour);
						} 
						// Add contour if it is arch shaped
						else if (contourIsArch(contour)) {
							prunedContours.push_back(contour);
						} 
						// Test if the contour can be transformed into an arch shape
						else {
							vector<Point> extractedArch = extractArch(contour);
							if (extractedArch != vector<Point>()) {
								prunedContours.push_back(extractedArch);
							}
						}
					} 
					// Add all contours
					else {
						prunedContours.push_back(contour);
					}
				}
			}
		}
	}
	return prunedContours;
}

vector<Point> extractArch(vector<Point> contour) {
	vector<Point> extractedArch = contour;
	Rect contourRec = boundingRect(contour);
	if (contourRec.width * 1.25f < contourRec.height) {
		// If 
		while (extractedArch.size() > 5) {
			Point lowestPoint = Point(0,0);
			int lowestIdx = -1;
			for (int i = 0; i < extractedArch.size();i++) {
				if (extractedArch[i].y > lowestPoint.y) {
					lowestPoint = extractedArch[i];
					lowestIdx = i;
				}
			}
			vector<Point> extractedArchTemp;
			for (int i = 0; i < extractedArch.size(); i++) {
				if (i != lowestIdx) {
					extractedArchTemp.push_back(extractedArch[i]);
				}
			}
			extractedArch = extractedArchTemp;
			if (contourIsArch(extractedArchTemp)) {
				return extractedArchTemp;
			}
		}
		if (contourIsArch(extractedArch)) {
			return extractedArch;
		}
		else {
			return vector<Point>();
		}
	}
	else {
		return vector<Point>();
	}
}

/** @brief Check if two points are close to each other 

@param p1 First Point
@param p2 Second Point
@param maxDeviation Maximal deviation in percent how much points can differ
*/
bool pointsAtSamePos(Point p1, Point p2, float maxDeviation) {
	int maxDistance = (min(innerFrameHeight, innerFrameWidth)) * maxDeviation;
	int distanceBetweenPoints = (int)norm(p1 - p2);
	if (distanceBetweenPoints < maxDistance) {
		return true;
	} else {
		return false;
	}
}

/*
extern "C" JNIEXPORT jstring JNICALL
Java_com_example_hellolibs_MainActivity_stringFromJNI(JNIEnv *env, jobject thiz) {
// Just for simplicity, we do this right away; correct way would do it in
// another thread...
#if defined(__arm__)
#if defined(__ARM_ARCH_7A__)
#if defined(__ARM_NEON__)
#if defined(__ARM_PCS_VFP)
#define ABI "armeabi-v7a/NEON (hard-float)"
#else
#define ABI "armeabi-v7a/NEON"
#endif
#else
#if defined(__ARM_PCS_VFP)
#define ABI "armeabi-v7a (hard-float)"
#else
#define ABI "armeabi-v7a"
#endif
#endif
#else
#define ABI "armeabi"
#endif
#elif defined(__i386__)
#define ABI "x86"
#elif defined(__x86_64__)
#define ABI "x86_64"
#elif defined(__mips64)  /* mips64el-* toolchain defines __mips__ too */
/*
#define ABI "mips64"
#elif defined(__mips__)
#define ABI "mips"
#elif defined(__aarch64__)
#define ABI "arm64-v8a"
#else
#define ABI "unknown"
#endif

return env->NewStringUTF("Hello from JNI !  Compiled with ABI " ABI ".");
}
*/
