#include <iostream>
#include <fstream>
#include <DRDSP/bitmap.h>

using namespace DRDSP;

Bitmap::Bitmap( uint32_t x, uint32_t y ) : sizeX(x), sizeY(y) {	
	pixels.resize( sizeX * sizeY );
	Clear();
}

void Bitmap::Clear() {
	memset(&pixels[0],0,sizeof(bmpPixel)*pixels.size());
}

void Bitmap::Clear( uint8_t R, uint8_t G, uint8_t B ) {
	bmpPixel rgb;
	rgb.b = B;
	rgb.g = G;
	rgb.r = R;
	for(auto& pixel : pixels) {
		pixel = rgb;
	}
}

void Bitmap::WriteFile( const char *filename ) const {
	ofstream out(filename,ios::out | ios::binary);
	if(!out) {
		cerr << "Bitmap::WriteFile -- Failed to open file " << filename << endl;
		return;
	}

	bmpFileHead fileHead( 54 + 3*sizeX*sizeY );
	bmpInfoHead infoHead( sizeX, sizeY );

	out.write( (const char*)&fileHead, sizeof(fileHead) );
	out.write( (const char*)&infoHead, sizeof(infoHead) );

	const char zeroBytes = 0x00;
	int numZeroBytes = (sizeX*sizeof(bmpPixel)) % 4;
	if( numZeroBytes > 0 )
		numZeroBytes = 4 - numZeroBytes;

	for(uint32_t j=0;j<sizeY;++j) {
		out.write( (const char*)&pixels[j*sizeX], sizeX*sizeof(bmpPixel) );
		out.write( &zeroBytes, numZeroBytes );
	}
}

void Bitmap::SetPixel( uint32_t x, uint32_t y, uint8_t R, uint8_t G, uint8_t B ) {
	bmpPixel& pixel = pixels[ x + (sizeY-y-1)*sizeX ];
	pixel.r = R;
	pixel.g = G;
	pixel.b = B;
}

const bmpPixel& Bitmap::operator()( uint32_t x, uint32_t y ) const {
	return pixels[ x + (sizeY-y-1)*sizeX ];
}

bmpPixel& Bitmap::operator()( uint32_t x, uint32_t y ) {
	return pixels[ x + (sizeY-y-1)*sizeX ];
}
