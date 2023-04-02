/*
 * Copyright (C) 2018, Ibon Guillen (http://giga.cps.unizar.es/~ibon/)
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom
 * the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _IMAGEIO_H_
#define _IMAGEIO_H_

#include "bunnykiller.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>
#include <string>
#include <vector>

/* Configure stb_image */

/* We make stb use the standard C++ allocators so we can
 pass around its buffers as-is */
namespace
{
	void* terrible_realloc(void* p, size_t oldsz, size_t newsz)
	{
		char* newp = new char[newsz]();
		std::copy((char*) p, (char*) p + oldsz, newp);

		delete[] ((char*) p); /* This is bad */

		return (void*) newp;
	}
	;
}

/* This is acceptable */
#define STBI_MALLOC(sz)                     new char[sz]()
/* This is bad */
#define STBI_FREE(p)                        delete[]((char*)p)
/* This is terrible */
#define STBI_REALLOC_SIZED(p, oldsz, newsz) terrible_realloc(p, oldsz, newsz)

#define STBI_NO_GIF /* The GIF loader requires an unsized realloc */

#define STB_IMAGE_IMPLEMENTATION
#include "External/stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "External/stb/stb_image_write.h"
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "External/stb/stb_image_resize.h"

#include "Image/Image.h"
#include "Image/RGBColor.h"
#include "Utils/Filesystem.h"

namespace Imaging
{
	template<class T>
	Image<T> load(const char* filename, unsigned load_channels = 0)
	{
		int w, h, c;
		std::unique_ptr<float[]> data(stbi_loadf(filename, &w, &h, &c, load_channels));

		if (!data) {
			throw std::runtime_error(
					"Error: failed to read image: \"" + std::string(filename) + "\"");
		}

		return Image<T>(std::move(data), w, h, c);
	}

	namespace
	{
		/* If the image is using internally a type different from float, we need to
		 * convert the data before storing it */
		template<class T>
		class Conversor
		{
		private:
			std::unique_ptr<float[]> m_data;
			size_t m_size;
		public:
			Conversor(const T* data, size_t size) :
					m_data(new float[size]()),
					m_size(size)
			{
				std::copy_n(data, size, m_data.get());
			}

			inline const float* data() const
			{
				return m_data.get();
			}

			inline size_t size() const
			{
				return m_size;
			}
		};

		/* If the image is using internally float, we can return it as-is */
		template<>
		class Conversor<float>
		{
		private:
			const float* m_data;
			size_t m_size;
		public:
			Conversor(const float* data, size_t size) :
					m_data(data),
					m_size(size)
			{
			}

			inline const float* data() const
			{
				return m_data;
			}

			inline size_t size() const
			{
				return m_size;
			}
		};
	}

	template<class T>
	void save(const Image<T>& img, const char* filename)
	{
		if (img.channels() > 3) {
			throw std::runtime_error(
					"Error: trying to save a " + std::to_string(img.channels())
							+ " channel image. Max supported output channels is 3");
		}

		/* Autodetect file format */
		std::string ext = Filesystem::file_extension(filename);

		int res = 1;
		if (ext == "hdr") {
			/* Save RGBE file */
			Conversor<T> conversor(img.data(), img.size());
			res = stbi_write_hdr(filename, (int) img.width(), (int) img.height(), (int) img.channels(), conversor.data());
		} else if (ext == "raw") {
			/* Save binary dump */
			FILE* file = std::fopen(filename, "w");
			size_t count = std::fwrite((const void*) img.data(), sizeof(T), img.size(),
					(FILE*) file);
			std::fclose(file);
			res = (count == img.size());
		} else if (ext == "png") {
			unsigned char *data = new unsigned char[(int) img.width() * (int) img.height() * 3];
			for(int i = 0; i < (int)img.height(); i++)
				for(int j = 0; j < (int)img.width(); j++)
				{
					data[i * img.width() * 3 + j * 3 + 0] = (unsigned char)(256 * img(j, i, 0));
					data[i * img.width() * 3 + j * 3 + 1] = (unsigned char)(256 * img(j, i, 1));
					data[i * img.width() * 3 + j * 3 + 2] = (unsigned char)(256 * img(j, i, 2));
				}
			res = stbi_write_png("abc.png", (int) img.width(),
					(int) img.height(), (int) img.channels(), data, (int) img.width() * 3);
			delete data;
		}
		else {
			throw std::runtime_error("Error: unsupported output file format: \"" + ext + "\"");
		}

		if (!res) {
			throw std::runtime_error(
					"Error: failed to write to output file: \"" + std::string(filename) + "\"");
		}
	}
}; /* namespace Imaging */

#endif /* _IMAGEIO_H_ */
