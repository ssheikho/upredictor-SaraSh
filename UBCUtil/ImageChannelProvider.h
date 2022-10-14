#ifndef IMAGE_CHANNEL_PROVIDER_H
#define IMAGE_CHANNEL_PROVIDER_H


#include <opencv2/opencv_modules.hpp>
#include <opencv2/opencv.hpp>
//#include <opencv/cv.hpp>

class ImageChannelProvider {
public:
	virtual cv::Mat &getImageChannel(size_t source, size_t index) = 0;
};

#endif
