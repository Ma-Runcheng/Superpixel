using TyImages
using Random
using Statistics
using TyPlot

function compute_densities_matrix(image, kernel_size)
    h, w, channels = size(image)
    inv_kernel_size_sqr = - 0.5 / (kernel_size * kernel_size)
    kernel_width = ceil(Int, 3 * kernel_size)
    densities = zeros(Float32, h, w)
    densities += randn(MersenneTwister(1234), Float32, (h, w)) * 0.00001
    for r = 1:h
        r_min = max(r - kernel_width, 1)
        r_max = min(r + kernel_width, h)
        for c = 1:w
            c_min = max(c - kernel_width, 1)
            c_max = min(c + kernel_width, w)
            for r_ = r_min:r_max
                for c_ = c_min:c_max
                    dist = 0.0
                    for channel = 1:channels
                        t = image[r, c, channel] - image[r_, c_, channel]
                        dist += t * t
                    end
                    t = r - r_
                    dist += t * t
                    t = c - c_
                    dist += t * t
                    densities[r, c] += exp(dist * inv_kernel_size_sqr)
                end
            end
        end
    end
    return densities
end


function compute_medoids_parent(image, dist_matrix, kernel_size)
    h, w, channels = size(image)
    kernel_width = ceil(Int, 3 * kernel_size)
    parent = reshape(collect(1:h*w), (h, w))  # 初始化parent为自己
    dist_parent = zeros(Float32, h, w)
    for r = 1:h
        r_min = max(r - kernel_width, 1)
        r_max = min(r + kernel_width, h)
        for c = 1:w
            current_density = dist_matrix[r, c]
            closest = typemax(Float32)
            c_min = max(c - kernel_width, 1)
            c_max = min(c + kernel_width, w)
            for r_ = r_min:r_max
                for c_ = c_min:c_max
                    if dist_matrix[r_, c_] > current_density
                        dist = 0.0
                        for channel = 1:channels
                            t = image[r, c, channel] - image[r_, c_, channel]
                            dist += t * t
                        end
                        t = r - r_
                        dist += t * t
                        t = c - c_
                        dist += t * t
                        if dist < closest
                            closest = dist
                            parent[r, c] = (c_ - 1) * w + r_
                        end
                    end
                end
            end
            dist_parent[r, c] = sqrt(closest)
        end
    end
    return reshape(parent, w * h), reshape(dist_parent, w * h)
end


function QuickShift(image, kernel_size, max_dist)
    h, w, channels = size(image)
    dens_matrix = compute_densities_matrix(image, kernel_size)  # (h, w)
    parent_flat, dist_parent_flat = compute_medoids_parent(image, dens_matrix, kernel_size)
    too_far = dist_parent_flat .> max_dist
    parent_flat[too_far] .= collect(1: h*w)[too_far]
    old = zeros(Int, h*w)

    while any(old .!= parent_flat)
        old .= parent_flat
        parent_flat .= parent_flat[parent_flat]
    end
    unique_index = unique(parent_flat)
    parent_flat = indexin(parent_flat, unique_index)
    parent = reshape(parent_flat, (h, w))
    return parent, dist_parent_flat
end

ratio = 0.2
origin_image = imread("image_path")
reshape_image = imresize(origin_image, (256, 256))
image = convert(Array{Float32}, reshape_image)
h = fspecial("gaussian",3,0)
image = imfilter(image, h)
image .= image * ratio
segmentation_mask, dist_parent_flat = QuickShift(image, 4, 200)

h = fspecial("laplacian");
segmentation_mask = segmentation_mask * 255 / maximum(segmentation_mask) 
boundary = imfilter(segmentation_mask, h)

for i = 1:size(image)[1]
    for j = 1:size(image)[2]
        if boundary[i, j] >= 1
            reshape_image[i, j, :] .= 1
        end
    end
end
imshow(reshape_image)
