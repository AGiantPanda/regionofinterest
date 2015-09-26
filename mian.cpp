#include "CImg.h"
//#include "contours.hpp"
#include "contour2.hpp"
#include <typeinfo>

using namespace std;
using namespace cimg_library;
using namespace roi;

#define IMGFILE "trachea_RGB.raw"
#define XMAX 370
#define YMAX 370
#define col 3
#define TYPE unsigned char

template<typename T> T *RGBtoPlanar(T *origin, int imgsize){
	T *tmp = new T[imgsize * 3];
	for (int i = 0; i < imgsize; i++){
		tmp[imgsize * 0 + i] = origin[i * 3 + 0];
		tmp[imgsize * 1 + i] = origin[i * 3 + 1];
		tmp[imgsize * 2 + i] = origin[i * 3 + 2];
	}
	return tmp;
}

template<typename T> 
T *FileRead(const char* FileName, int width, int height, int col_t = 1)
{
	FILE *pFile = std::fopen(FileName, "rb");
	if (NULL == pFile){
		std::cout << "read failed1" << std::endl;
		exit(0);
	}
	int imgsize = width*height*col_t;
	T *pImg = new T[imgsize];
	short *type;
	unsigned short *type2;
	bool readok = (imgsize == fread(pImg, sizeof(T), imgsize, pFile));
	if (!readok){
		std::cout << "read failed2" << std::endl;
		exit(0);
	}
	std::fclose(pFile);

	//大小端转换以及强度最大最小值映射
	if (typeid(pImg) == typeid(type) || typeid(pImg) == typeid(type2)) {
		int low = 1000, high = 0;
		for (int i = 0; i < imgsize; i++){
			pImg[i] = (pImg[i] << 8) | ((pImg[i] >> 8) & 0xFF);
			low = pImg[i] < low ? pImg[i] : low;
			high = pImg[i]>high ? pImg[i] : high;
		}
		for (int i = 0; i < imgsize; i++){
			pImg[i] = (pImg[i] - low) / ((high - low) / 255);
		}
		std::cout << low << std::endl << high << std::endl;
	}
	//else if (col_t == 3){
	//	T *pImg3 = new T[imgsize];
	//	for (int i = 0; i < XMAX*YMAX; i++){
	//		pImg3[XMAX*YMAX * 0 + i] = pImg[i * 3 + 0];
	//		pImg3[XMAX*YMAX * 1 + i] = pImg[i * 3 + 1];
	//		pImg3[XMAX*YMAX * 2 + i] = pImg[i * 3 + 2];
	//	}
	//	return pImg3;
	//}

	return pImg;
}

int main(int argc, char* argv[])
{
	int imgsize = XMAX * YMAX;
	TYPE* pImg = FileRead<TYPE>(IMGFILE, XMAX, YMAX, col);

	CImg<TYPE> raw;
	if (col == 3){
		CImg<unsigned char> tmp(RGBtoPlanar<TYPE>(pImg, XMAX*YMAX), XMAX, YMAX, 1, col);
		raw = tmp;
	}
	else{
		CImg<unsigned char> tmp(pImg, XMAX, YMAX, 1, col);
		raw = tmp;
	}
	CImgDisplay disp(raw, "rawpic");
	CImgDisplay disp1(XMAX, YMAX, "contour_pic", 0);
	CImgDisplay disp2(XMAX, YMAX, "original_pic", 0);
	CImgDisplay disp3(XMAX, YMAX, "scanline_pic", 0);
	disp2.move(400, 100);
	disp3.move(800, 100);
	//CImgDisplay bi_disp(512, 512, "bi_pic", 0);
	int mx = 0, my = 0; 
	bool redraw = false;
	bool move = false;
	vector<Point> points;
	const unsigned char white[] = { 255, 255, 255 };
	(+raw).draw_text(1, 1, "press right button to start.", white, 0, 1, 20).display(disp);

	while (!disp.is_closed() && !disp2.is_closed()){
		if (disp.mouse_x() >= 0){
			mx = disp.mouse_x(); my = disp.mouse_y();
		}
		if (disp.button() & 2){
			if (!move){
				points.clear();
				(+raw).draw_text(1, 1, "move & click, press right button to stop.\ndo not draw a circle with any cross point.", white, 0, 1, 20).display(disp);
				disp.set_key(); move = true;
			}
			else{
				(+raw).draw_text(1, 1, "press right button to start.", white, 0, 1, 20).display(disp);
				disp.set_key(); move = false; redraw = true;
			}
		}
		if ((disp.button() & 1) && move){
			points.push_back(Point(mx, my));
			if (points.size() > 1){
				CImg<TYPE> raw_tmp = raw;
				draw_line_loop<TYPE>(raw_tmp, points);
				(+raw_tmp).draw_text(1, 1, "move & click, press right button to stop.\ndo not draw a circle with any cross point.", white, 0, 1, 20).display(disp);
			}
			//disp.set_button();
		}
		if (redraw){
			//x0 = mx - factor; y0 = my - factor;
			//x1 = mx + factor; y1 = my + factor;
			//if (x0 < 0){
			//	x0 = 0; x1 = x0 + 2 * factor;
			//}
			//if (y0 < 0){
			//	y0 = 0; y1 = y0 + 2 * factor;
			//}
			//if (x1 > XMAX){
			//	x1 = XMAX; x0 = x1 - 2 * factor;
			//}
			//if (y1 > YMAX){
			//	y1 = YMAX; y0 = y1 - 2 * factor;
			//}
			//const unsigned char white[] = { 255, 255, 255 };
			////(+raw).draw_rectangle(x0, y0, x1, y1, white, 1.0f, ~0U).draw_text(2, 2, "x=%d, y=%d", white, 0, 1, 13, mx, my).draw_line(100, 100, mx, my, white).display(disp);
			//points.push_back(Point(x0, y0));
			//points.push_back(Point(x0, y1)); 
			//points.push_back(Point(x0 + 10 , y1));
			//points.push_back(Point(x0 + 25, y1 - 35));
			//points.push_back(Point(x0 + 40, y1 - 10));
			//points.push_back(Point(x0 + 60, y1 - 30));
			//points.push_back(Point(x0 + 80, y1 - 20));
			//points.push_back(Point(x1, y1));
			//points.push_back(Point(x1, y0));
			//points.push_back(Point(x0 + 80, y0));
			//points.push_back(Point(x0 + 65, y0 + 28));
			//points.push_back(Point(x0 + 90, y0 + 30));
			//points.push_back(Point(x0 + 60, y0 + 30));
			//points.push_back(Point(x0 + 40, y0));
			//points.push_back(Point(x0 + 25, y0 + 25));
			//points.push_back(Point(x0 + 10, y0));
			
			//CImg<unsigned char> bimask = AreaMask<unsigned char>(raw, points);
			//bimask.display(bi_disp);
			//cout << "+++++++++++++++++++++++++++++++++++" << endl;
			//CImg<unsigned char> bimask;
			//int bitime = clock();
			//int threshold = otsu(raw, points, bimask);
			//cout << "this otsu takes: " << clock() - bitime << endl;
			//bitime = clock();
			////the following code takes too much time
			//unsigned char *bi_pic = new unsigned char[raw.width()*raw.height()];
			//int height = raw.height(), width = raw.width();
			//for (int i = 0; i < height; i++){
			//	for (int j = 0; j < width; j++){
			//		if (raw.atXY(j, i) >= threshold && bimask.atXY(j, i) == 255)
			//			bi_pic[i*width + j] = 255;
			//		else
			//			bi_pic[i*width + j] = 0;
			//	}
			//}
			////返回二值图像
			//cimg_library::CImg<unsigned char> bis(bi_pic, raw.width(), raw.height());
			//delete[]bi_pic;
			//bis.display(bi_disp);
			//cout << "bi compute time is: " << clock() - bitime << endl << "---------------------------------" << endl;
			int begin = clock();
			int total = clock();
			vector<Point> pr = FindBiggestContour<TYPE>(raw, points, 1, 1);
			cout << "- number of points: " << points.size() << endl;
			cout << "- findbiggestcontour time: " << clock() - begin << endl;
			begin = clock();
			vector<Point>::iterator it;
			memset(pImg, 0, imgsize*sizeof(unsigned char));
			for (it = pr.begin(); it != pr.end(); it++)
				pImg[(*it).y * XMAX + (*it).x] = 255;

			CImg<unsigned char> bil(pImg, XMAX, YMAX);
			cout << "- generate pic time: " << clock() - begin << endl
				<< "TOTAL TIME: " << clock() - total << endl;
			(+bil).draw_text(2, 2, "total time: %d", white, 0, 1, 13, clock() - total).display(disp1);
			//int mask_time = clock();
			//bil = AreaMask(raw, points);
			//(+bil).draw_text(2, 2, "original time: %d", white, 0, 1, 13, clock() - mask_time).display(disp2);
			//mask_time = clock();
			//bil = _AreaMask(raw, points);
			//(+bil).draw_text(2, 2, "scanline time: %d", white, 0, 1, 13, clock() - mask_time).display(disp3);
			redraw = false;
		}
		CImgDisplay::wait(disp);
	}
	//vector<MYPOINT> pointresault = FindBiggestContour(raw, points);
	//vector<MYPOINT>::iterator it;
	//memset(pImg, 0, imgsize*sizeof(unsigned char));
	//for (it = pointresault.begin(); it != pointresault.end(); it++)
	//	pImg[(*it).x * 512 + (*it).y] = 255;
	//CImg<unsigned char> bil(pImg, 512, 512);
	//CImgDisplay disp2(bil, "dfa");

	//getchar(); getchar();
	delete[]pImg;
	return 0;
}