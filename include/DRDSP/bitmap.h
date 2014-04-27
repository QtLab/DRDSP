#ifndef INCLUDED_BITMAP
#define INCLUDED_BITMAP
#include <stdint.h>
#include <vector>

using namespace std;

namespace DRDSP {

	#pragma pack(1)
	struct bmpPixel {
		uint8_t b,g,r;
	};

	struct bmpFileHead {
		uint16_t bfType;
		uint32_t bfSize;
		uint16_t bfReserved1;
		uint16_t bfReserved2;
		uint32_t bfOffBits;

		bmpFileHead() :
			bfType(0x4D42),
			bfReserved1(0),
			bfReserved2(0),
			bfOffBits(54)
		{}

	};

	struct bmpInfoHead {
		uint32_t biSize;
		uint32_t biWidth;
		uint32_t biHeight;
		uint16_t biPlanes;
		uint16_t biBitCount;
		uint32_t biCompression;
		uint32_t biSizeImage;
		uint32_t biXPelsPerMeter;
		uint32_t biYPelsPerMeter;
		uint32_t biClrUsed;
		uint32_t biClrImportant;

		bmpInfoHead() :
			biSize(40),
			biPlanes(1),
			biBitCount(24),
			biCompression(0),
			biSizeImage(0),
			biXPelsPerMeter(0),
			biYPelsPerMeter(0),
			biClrUsed(0),
			biClrImportant(0)
		{}

	};
	#pragma pack()

	struct Bitmap {

		Bitmap(uint32_t x,uint32_t y);
		void Clear();
		void Clear( uint8_t R, uint8_t G, uint8_t B );
		void SetPixel( uint32_t x, uint32_t y, uint8_t R, uint8_t G, uint8_t B );
		const bmpPixel& operator()( uint32_t x, uint32_t y ) const;
		bmpPixel& operator()( uint32_t x, uint32_t y );
		void WriteFile( const char *filename ) const;
	
	protected:
		uint32_t sizeX, sizeY;
		vector<bmpPixel> pixels;
		bmpFileHead fileHeader;
		bmpInfoHead infoHeader;
	};

}

#endif
