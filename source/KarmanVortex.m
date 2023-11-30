clear; clc;

% --------------------------------------
% 各種定数
% --------------------------------------
height = 80;
width = 200;
viscosity = 0.02;
u0 = 0.10;
frame = 1000;
omega = 1 / (3*viscosity + 0.5);

% --------------------------------------
% objectの作成
% --------------------------------------
object = false(height, width);
object(floor(height/2)-8:floor(height/2)+8, floor(height/2)) = true;

% 乱れのために微小の非対称性を与える
object(floor(height/2)-8, floor(height/2) + 1) = true;

object_N  = circshift(object,    1, 1);
object_S  = circshift(object,   -1, 1);
object_E  = circshift(object,    1, 2);
object_W  = circshift(object,   -1, 2);
object_NE = circshift(object_N,  1, 2);
object_NW = circshift(object_N, -1, 2);
object_SE = circshift(object_S,  1, 2);
object_SW = circshift(object_S, -1, 2);

% --------------------------------------
% f_iの作成
% --------------------------------------
f_0  = 4.0/9.0  * ( ones(height,width)                   - 1.5*u0^2 );
f_N  = 1.0/9.0  * ( ones(height,width)                   - 1.5*u0^2 );
f_S  = 1.0/9.0  * ( ones(height,width)                   - 1.5*u0^2 );
f_E  = 1.0/9.0  * ( ones(height,width) + 3*u0 + 4.5*u0^2 - 1.5*u0^2 );
f_W  = 1.0/9.0  * ( ones(height,width) - 3*u0 + 4.5*u0^2 - 1.5*u0^2 );
f_NE = 1.0/36.0 * ( ones(height,width) + 3*u0 + 4.5*u0^2 - 1.5*u0^2 );
f_SE = 1.0/36.0 * ( ones(height,width) + 3*u0 + 4.5*u0^2 - 1.5*u0^2 );
f_NW = 1.0/36.0 * ( ones(height,width) - 3*u0 + 4.5*u0^2 - 1.5*u0^2 );
f_SW = 1.0/36.0 * ( ones(height,width) - 3*u0 + 4.5*u0^2 - 1.5*u0^2 );

rho = f_0 + f_N + f_S + f_E + f_W + f_NE + f_SE + f_NW + f_SW;
u = (f_E + f_NE + f_SE - f_W - f_NW - f_SW) ./ rho;
v = (f_N + f_NE + f_NW - f_S - f_SE - f_SW) ./ rho;


% --------------------------------------
% 実行部分
%   - stream, collide
%   - 描画対象の作成
%   - 描画
% --------------------------------------

frames(frame) = struct('cdata', [], 'colormap', []);
fig = figure;
for i = 1:frame
    % stream, collide
    for j = 1:20
        [object, f_N, f_S, f_E, f_W, f_NE, f_NW, f_SE, f_SW, object_N, object_S, object_E, object_W, object_NE, object_NW, object_SE, object_SW] = stream(object, f_N, f_S, f_E, f_W, f_NE, f_NW, f_SE, f_SW, object_N, object_S, object_E, object_W, object_NE, object_NW, object_SE, object_SW);
        [rho, u0, u, v, f_0, f_N, f_S, f_E, f_W, f_NE, f_NW, f_SE, f_SW, omega] = collide(rho, u0, u, v, f_0, f_N, f_S, f_E, f_W, f_NE, f_NW, f_SE, f_SW, omega);
    end
    
    % 描画対象の作成
    curl_result = curl(u, v);

    maxim = max(curl_result,[],'all');
    minim = min(curl_result,[],'all');

    curl_result = (curl_result / maxim + 0.5) * 255;

    ObjectArray = 255 * object;
    curl_result = curl_result + ObjectArray;
    
    % 描画
    clf;
    image(curl_result);
    colormap('jet')
    pbaspect([width height 1])
    xlim([0 width]); ylim([0 height]);
    title(['step : ', num2str(i)])
    drawnow;
    frames(i) = getframe(fig);
end

%% output
exporttype = 'mp4';
% exporttype = 'gif';
switch exporttype
    case 'mp4'
        video = VideoWriter('karmanvortex.mp4', 'MPEG-4');
        open(video);
        writeVideo(video, frames);
        close(video);
    case 'gif'
        filename = 'KarmanVortex.gif';
        for i = 1:frame
            [A, map] = rgb2ind(frame2im(frames(i)), 256);
            if i == 1
                imwrite(A, map, filename, 'gif', 'DelayTime', 1/30);
            else
                imwrite(A, map, filename, 'gif', 'DelayTime', 1/30, 'WriteMode', 'append');
            end
        end
end

%% function
% --------------------------------------
%   - stream
%   - collide
%   - curl
% --------------------------------------

function [rho, u0, u, v, f_0, f_N, f_S, f_E, f_W, f_NE, f_NW, f_SE, f_SW, omega] = collide(~, u0, ~, ~, f_0, f_N, f_S, f_E, f_W, f_NE, f_NW, f_SE, f_SW, omega)
    rho = f_0 + f_N + f_S + f_E + f_W + f_NE + f_SE + f_NW + f_SW;
    u = (f_E + f_NE + f_SE - f_W - f_NW - f_SW) ./ rho;
    v = (f_N + f_NE + f_NW - f_S - f_SE - f_SW) ./ rho;
    ux2 = u .* u;
    uy2 = v .* v;
    uxuy = u .* v;
    u2 = ux2 + uy2;
    
    f_0  = (1 - omega) * f_0  + omega *  4.0/9.0 * rho .* ( 1                                        - 1.5 * u2);
    f_N  = (1 - omega) * f_N  + omega *  1.0/9.0 * rho .* ( 1 + 3 * (     v) + 4.5 * (          uy2) - 1.5 * u2);
    f_S  = (1 - omega) * f_S  + omega *  1.0/9.0 * rho .* ( 1 + 3 * (    -v) + 4.5 * (          uy2) - 1.5 * u2);
    f_E  = (1 - omega) * f_E  + omega *  1.0/9.0 * rho .* ( 1 + 3 * (     u) + 4.5 * (          ux2) - 1.5 * u2);
    f_W  = (1 - omega) * f_W  + omega *  1.0/9.0 * rho .* ( 1 + 3 * (    -u) + 4.5 * (          ux2) - 1.5 * u2);
    f_NE = (1 - omega) * f_NE + omega * 1.0/36.0 * rho .* ( 1 + 3 * ( u + v) + 4.5 * (u2 + 2 * uxuy) - 1.5 * u2);
    f_NW = (1 - omega) * f_NW + omega * 1.0/36.0 * rho .* ( 1 + 3 * (-u + v) + 4.5 * (u2 - 2 * uxuy) - 1.5 * u2);
    f_SE = (1 - omega) * f_SE + omega * 1.0/36.0 * rho .* ( 1 + 3 * ( u - v) + 4.5 * (u2 - 2 * uxuy) - 1.5 * u2);
    f_SW = (1 - omega) * f_SW + omega * 1.0/36.0 * rho .* ( 1 + 3 * (-u - v) + 4.5 * (u2 + 2 * uxuy) - 1.5 * u2);

    f_E(:,1)   = 1.0/9.0  * (1 + 3 *  u0 + 4.5 * u0^2 - 1.5 * u0^2);
    f_W(:, 1)  = 1.0/9.0  * (1 + 3 * -u0 + 4.5 * u0^2 - 1.5 * u0^2);
    f_NE(:, 1) = 1.0/36.0 * (1 + 3 *  u0 + 4.5 * u0^2 - 1.5 * u0^2);
    f_SE(:, 1) = 1.0/36.0 * (1 + 3 *  u0 + 4.5 * u0^2 - 1.5 * u0^2);
    f_NW(:, 1) = 1.0/36.0 * (1 + 3 * -u0 + 4.5 * u0^2 - 1.5 * u0^2);
    f_SW(:, 1) = 1.0/36.0 * (1 + 3 * -u0 + 4.5 * u0^2 - 1.5 * u0^2);
end

function [object, f_N, f_S, f_E, f_W, f_NE, f_NW, f_SE, f_SW, object_N, object_S, object_E, object_W, object_NE, object_NW, object_SE, object_SW] = stream(object, f_N, f_S, f_E, f_W, f_NE, f_NW, f_SE, f_SW, object_N, object_S, object_E, object_W, object_NE, object_NW, object_SE, object_SW)
    f_N = circshift(f_N,  1, 1);
    f_S = circshift(f_S, -1, 1);
    f_E = circshift(f_E,  1, 2);
    f_W = circshift(f_W, -1, 2);

    f_NE = circshift(f_NE,  1, 1);
    f_NE = circshift(f_NE,  1, 2);
    f_NW = circshift(f_NW,  1, 1);
    f_NW = circshift(f_NW, -1, 2);
    f_SE = circshift(f_SE, -1, 1);
    f_SE = circshift(f_SE,  1, 2);
    f_SW = circshift(f_SW, -1, 1);
    f_SW = circshift(f_SW, -1, 2);

    f_N(object_N)   = f_S(object);
    f_S(object_S)   = f_N(object);
    f_E(object_E)   = f_W(object);
    f_W(object_W)   = f_E(object);
    f_NE(object_NE) = f_SW(object);
    f_NW(object_NW) = f_SE(object);
    f_SE(object_SE) = f_NW(object);
    f_SW(object_SW) = f_NE(object);
end

function curl_result = curl(u, v)
    curl_result = circshift(u, 1, 1) - circshift(u, -1, 1) + circshift(v, -1, 2) - circshift(v, 1, 2);
end