#include <iostream>
#include <fstream>
#include <DRDSP/bitmap.h>

using namespace DRDSP;

Bitmap::Bitmap( uint32_t x, uint32_t y ) : sizeX(x), sizeY(y) {	
	pixels.resize( sizeX * sizeY );

	fileHeader.bfSize = 54 + 3*sizeX*sizeY;
	infoHeader.biWidth = sizeX;
	infoHeader.biHeight = sizeY;

	Clear();
}

void Bitmap::Clear() {
	memset(&pixels[0],0,sizeof(bmpPixel)*pixels.size());
}

void Bitmap::Clear( uint8_t R, uint8_t G, uint8_t B ) {
	for(auto& pixel : pixels) {
		pixel.b = B;
		pixel.g = G;
		pixel.r = R;
	}
}

void Bitmap::WriteFile( const char *filename ) const {
	ofstream out(filename,ios::out | ios::binary); //Access outfile

	if(!out) {
		cerr << "Bitmap::WriteFile -- Failed to open file " << filename << endl;
		return;
	}

	out.write((const char*) &fileHeader,sizeof(fileHeader));
	out.write((const char*) &infoHeader,sizeof(infoHeader));

	char zeroBytes = 0x00;
	uint8_t numZeroBytes = ( (sizeX*sizeof(bmpPixel)) % 4 );
	if( numZeroBytes > 0 )
		numZeroBytes = 4 - numZeroBytes;

	for(uint32_t j=0;j<sizeY;++j) {
		out.write( (const char*)&pixels[j*sizeX], sizeX*sizeof(bmpPixel) );
		out.write( (const char*)&zeroBytes, numZeroBytes );
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
