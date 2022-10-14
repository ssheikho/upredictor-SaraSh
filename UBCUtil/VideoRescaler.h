#ifndef VIDEO_RESCALER_H
#define VIDEO_RESCALER_H

#include "TypeBroadcaster.h"
#include "TypeRecipient.h"


//#include <opencv/cv.hpp>


#include <opencv2/opencv_modules.hpp>
#include <opencv2/opencv.hpp>
#include "opencv2/core/core_c.h"
#include "opencv2/videoio/legacy/constants_c.h"
#include "opencv2/highgui/highgui_c.h"


class VideoRescaler :
	public TypeBroadcaster<cv::Mat>
	, public TypeRecipient<cv::Mat> {
public:
	VideoRescaler(int destW, int destH);
	~VideoRescaler();

	void processType(const cv::Mat t);

protected:
	cv::Mat _out;
};

#endif
