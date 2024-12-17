using TyImages
using Random
using Statistics
using TyPlot
using ImageView, ImageCore, BenchmarkTools
using Images

mutable struct Cluster
    l
    a
    b
    y
    x
end

function init_centroids(img_lab, S)
    h, w = size(img_lab)
    clusters = Cluster[]
    labels = fill(-1, h, w)
    pixel_count = Integer[]
    for x = div(S, 2):S:w
        for y=div(S, 2):S:h
            push!(clusters, Cluster(img_lab[y, x].l,
                                   img_lab[y, x].a,
                                   img_lab[y, x].b,
                                   y,
                                   x))
            push!(pixel_count, 0)
        end
    end
    return clusters, labels, pixel_count
end

function moveCenter(clusters, img_lab)
    h, w = size(img_lab)
    function get_gradient(y, x)
        if x + 1 > w x = w - 2 end
        if y + 1 > h y = h - 2 end
    
        return img_lab[y + 1, x + 1].l - img_lab[y, x].l + 
               img_lab[y + 1, x + 1].a - img_lab[y, x].a + 
               img_lab[y + 1, x + 1].b - img_lab[y, x].b
    end

    for i = 1:length(clusters)
        current_grad = get_gradient(clusters[i].y, clusters[i].x)
        for dh = -1:1
            for dw = -1:1
                _y = clusters[i].y + dh
                _x = clusters[i].x + dw
                new_gradient = get_gradient(_y, _x)
                if new_gradient < current_grad
                    clusters[i].l = img_lab[_y, _x].l
                    clusters[i].a = img_lab[_y, _x].a
                    clusters[i].b = img_lab[_y, _x].b
                    clusters[i].y = _y
                    clusters[i].x = _x
                    current_grad = new_gradient 
                end
            end
        end
    end
    return clusters
end

function cluters_pixels(clusters, labels, S, M,img_lab)
    h, w = size(img_lab)
    distance = fill(Inf, h, w)  # 初始化距离种子点距离无穷大
    for i = 1:length(clusters)  # 遍历2S*2S范围内的邻域
        for x = (clusters[i].x - 2 * S):(clusters[i].x + 2 * S)
            if x <= 0 || x > w continue end
            for y = (clusters[i].y - 2 * S):(clusters[i].y + 2 * S)
                if y <= 0 || y > h continue end

                L = img_lab[y, x].l
                A = img_lab[y, x].a
                B = img_lab[y, x].b
                Dc = sqrt((L - clusters[i].l)^2 + 
                          (A - clusters[i].a)^2 +
                          (B - clusters[i].b)^2)
                Ds = sqrt((y - clusters[i].y)^2 +
                          (x - clusters[i].x)^2)
                D = sqrt((Dc / M)^2 + (Ds / S)^2)

                if D < distance[y, x]
                    distance[y, x] = D
                    labels[y, x] = i
                end
            end
        end
    end
    return clusters, labels
end

function update_cluster_position(clusters, labels, pixel_count, img_lab)
    h, w = size(img_lab)
    for i = 1:length(clusters)
        clusters[i].y = clusters[i].x = pixel_count[i] = 0  # 将所有种子点位置设为0
    end
    for x in 1:w
        for y in 1:h
            label_index = labels[y, x]
            if label_index == -1 continue end

            clusters[label_index].y += y
            clusters[label_index].x += x
            pixel_count[label_index] += 1
        end
    end
    # 取聚类区域内的均值
    for i = 1:length(clusters)
        new_y = div(clusters[i].y, pixel_count[i])
        new_x = div(clusters[i].x, pixel_count[i])
        clusters[i].l = img_lab[new_y, new_x].l
        clusters[i].a = img_lab[new_y, new_x].a
        clusters[i].b = img_lab[new_y, new_x].b
        clusters[i].y = new_y
        clusters[i].x = new_x
    end

    return clusters, labels
end

function SLIC(image, n_segments,M, step)
    image_Lab = Lab.(image)  # 将图片转换为Lab色彩空间
    h, w = size(image)
    S = round(Int, (sqrt((h * w) / n_segments)))  # S为间隔
    clusters, labels, pixel_count = init_centroids(image_Lab, S)  # 初始化中心
    clusters = moveCenter(clusters, image_Lab)  # 重选种子点
    # 迭代优化
    for i = 1:step
        print(i)
        clusters, labels = cluters_pixels(clusters, labels, S, M, image_Lab)
        clusters, labels = update_cluster_position(clusters, labels, pixel_count, image_Lab)
    end
    return labels
end

image = load("image_path")

step = 10
n_segments = 100
M = 100
labels = SLIC(image, n_segments, M, step)

h = fspecial("laplacian");
segmentation_mask = labels * 255 / maximum(labels) 
boundary = TyImages.imfilter(segmentation_mask, h)
for i = 1:size(image)[1]
    for j = 1:size(image)[2]
        if boundary[i, j] >= 1
            image[i, j, :] .= 1
        end
    end
end
TyImages.imshow(image)