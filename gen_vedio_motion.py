import cv2
import numpy as np
import glob
import matplotlib.pyplot as plt

# 读取HDR图像
image_files = sorted(glob.glob('./build/Release/cornell_*.hdr'))
# image_files = sorted(glob.glob('./build/Release/cube_*.hdr'))

images = []
for file in image_files:
    image = cv2.imread(file, cv2.IMREAD_ANYDEPTH | cv2.IMREAD_COLOR)
    # plt.imshow(image)
    # plt.show()
    images.append(image)

images = np.array(images) # h,w,t,c
images[images > 0.999] = 0.999
print(images.shape)

# 创建视频
output_path = 'output.mp4'
fourcc = cv2.VideoWriter_fourcc(*"mp4v")
video_writer = cv2.VideoWriter(output_path, fourcc, 25.0, (images.shape[2], images.shape[1]))

# 写入视频
images = np.power(images, 1./4.)
for image in images:
    image = (image * 256).astype(np.int8)
    video_writer.write(image)

# 释放资源
video_writer.release()
