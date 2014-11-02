#ifndef INCLUDED_BITMAP
#define INCLUDED_BITMAP
#include <cstdint>
#include <vector>

using namespace std;

namespace DRDSP {

	#pragma pack(1)
	struct bmpPixel {
		uint8_t b,g,r;
	};

	struct bmpFileHead {
		uint16_t bfType = 0x4D42;
		uint32_t bfSize;
		uint16_t bfReserved1 = 0;
		uint16_t bfReserved2 = 0;
		uint32_t bfOffBits = 54;

		bmpFileHead() = default;
		explicit bmpFileHead( uint32_t size ) : bfSize(size) {}
	};

	struct bmpInfoHead {
		uint32_t biSize = 40;
		uint32_t biWidth;
		uint32_t biHeight;
		uint16_t biPlanes = 1;
		uint16_t biBitCount = 24;
		uint32_t biCompression = 0;
		uint32_t biSizeImage = 0;
		uint32_t biXPelsPerMeter = 0;
		uint32_t biYPelsPerMeter = 0;
		uint32_t biClrUsed = 0;
		uint32_t biClrImportant = 0;

		bmpInfoHead() = default;
		bmpInfoHead( uint32_t width, uint32_t height ) : biWidth(width), biHeight(height) {}
	};
	#pragma pack()

	struct Bitmap {
		Bitmap(uint32_t x,uint32_t y);
		void Clear();
		void Clear( uint8_t R, uint8_t G, uint8_t B );
		void SetPixel( uint32_t x, uint32_t y, uint8_t R, uint8_t G, uint8_t B );
		const bmpPixel& operator()( uint32_t x, uint32_t y ) const;
		bmpPixel& operator()( uint32_t x, uint32_t y );
		void WriteFile( const char* filename ) const;
	
	protected:
		uint32_t sizeX, sizeY;
		vector<bmpPixel> pixels;
	};

}

#endif
