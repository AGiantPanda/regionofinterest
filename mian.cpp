#include "CImg.h"
#include "contours.hpp"

using namespace std;
using namespace cimg_library;

#define XMAX 512
#define YMAX 512

int main(int argc, char* argv[]){
	FILE *pFile = fopen("512_512.raw", "rb");
	int imgsize = 512 * 512;
	unsigned char* pImg = new unsigned char[imgsize];
	bool readof = (imgsize == fread(pImg, sizeof(unsigned char), imgsize, pFile));
	fclose(pFile);

	CImg<unsigned char> raw(pImg, 512, 512);
	CImgDisplay disp(raw, "rawpic");
	CImgDisplay disp2(512, 512, "dd", 0);
	int mx = 0, my = 0, factor = 50, x0, y0, x1, y1; 
	bool redraw = false;

	while (!disp.is_closed() && !disp2.is_closed()){
		if (disp.mouse_x() >= 0){
			mx = disp.mouse_x(); my = disp.mouse_y(); redraw = true;
		}
		if (redraw){
			x0 = mx - factor; y0 = my - factor;
			x1 = mx + factor; y1 = my + factor;
			if (x0 < 0){
				x0 = 0; x1 = x0 + 2 * factor;
			}
			if (y0 < 0){
				y0 = 0; y1 = y0 + 2 * factor;
			}
			if (x1 > XMAX){
				x1 = XMAX; x0 = x1 - 2 * factor;
			}
			if (y1 > YMAX){
				y1 = YMAX; y0 = y1 - 2 * factor;
			}
			const unsigned char white[] = { 255, 255, 255 };
			(+raw).draw_rectangle(x0, y0, x1, y1, white, 1.0f, ~0U).display(disp);

			vector<MYPOINT> points;
			points.push_back(MYPOINT(x0, y0)); points.push_back(MYPOINT(x1, y0)); points.push_back(MYPOINT(x1, y1)); points.push_back(MYPOINT(x0, y1));
			vector<MYPOINT> pr = FindBiggestContour(raw, points);
			vector<MYPOINT>::iterator it;
			memset(pImg, 0, imgsize*sizeof(unsigned char));
			for (it = pr.begin(); it != pr.end(); it++)
				pImg[(*it).y * 512 + (*it).x] = 255;
			/*int th = otsu(raw, points, 1);
			for (int i = 0; i < raw.height();i++)
			for (int j = 0; j < raw.width(); j++){
				if (raw.atXY(j, i) >= th) pImg[i*raw.width() + j] = 255;
				else pImg[i*raw.width() + j] = 0;
			}*/
			CImg<unsigned char> bil(pImg, 512, 512);
			bil.display(disp2);
		}
		if (disp.button() & 1){
			factor = (int)(factor / 1.5f);
			if (factor<3) factor = 3;
			disp.set_button(); redraw = true;
		}
		if (disp.button() & 2){
			factor = (int)(factor * 1.5f);
			if (factor>100) factor = 100;
			disp.set_button(); redraw = true;
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

	getchar(); getchar();
	return 0;
}