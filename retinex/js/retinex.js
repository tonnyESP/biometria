//Retinex
//(c) 2013 by Tomasz Lubinski
//www.algorytm.org
/* Data of the image */
var imageData;

/* Papers:  "Recursive Implementation of the gaussian filter.",
Ian T. Young , Lucas J. Van Vliet, Signal Processing 44, Elsevier 1995. */

function GetGaussianCoeficients(sigma) {
    var result = new Array();
    var q = 0;
    if (sigma >= 2.5) q = 0.98711 * sigma - 0.96330;
    else if ((sigma >= 0.5) && (sigma < 2.5)) q = 3.97156 - 4.14554 * Math.sqrt(1.0 - 0.26891 * sigma);
    else
        q = 0.1147705018520355224609375;

    var q2 = q * q;
    var q3 = q * q2;
    result.push(1.57825 + (2.44413 * q) + (1.4281 * q2) + (0.422205 * q3));
    result.push((2.44413 * q) + (2.85619 * q2) + (1.26661 * q3));
    result.push(-((1.4281 * q2) + (1.26661 * q3)));
    result.push((0.422205 * q3));
    result.push(1.0 - ((result[1] + result[2] + result[3]) / result[0]));
    return result;
}

function GaussianSmooth(imgIn, imgOut, shift, rowstride, size, coefs) {
    // forward pass
    var w1 = new Array();
    var w2 = new Array();
    w1[0] = imgIn[shift];
    w1[1] = imgIn[shift];
    w1[2] = imgIn[shift];
    for (var i = 0, n = 3; i < size; i++, n++) {
        w1[n] = (coefs[4] * imgIn[shift + i * rowstride] + ((coefs[1] * w1[n - 1] + coefs[2] * w1[n - 2] + coefs[3] * w1[n - 3]) / coefs[0]));
    }
    w1[size + 0] = w1[size - 1];
    w1[size + 1] = w1[size - 2];
    w1[size + 2] = w1[size - 3];

    // backward pass
    w2[size + 0] = w1[size + 0];
    w2[size + 1] = w1[size + 1];
    w2[size + 2] = w1[size + 2];
    w2[size + 3] = w1[size + 2];
    w2[size + 4] = w1[size + 2];
    w2[size + 5] = w1[size + 2];
    for (var i = size - 1, n = i + 3; i >= 0; i--, n--) {
        w2[n] = imgOut[shift + i * rowstride] = (coefs[4] * w1[n] + ((coefs[1] * w2[n + 1] + coefs[2] * w2[n + 2] + coefs[3] * w2[n + 3]) / coefs[0]));
    }
}

/* Retinex */
function Retinex(sigmas, weights, alfa, beta) 
{
    var canvas = document.getElementById("myCanvas");
    var ctx = canvas.getContext("2d");

    //read the width and height of the canvas
    var width = canvas.width;
    var height = canvas.height;

    var bluredImageData = ctx.createImageData(width, height).data;
    var imageDataData = imageData.data;

    //initialize result
    var retinexResult = new Array();
    //for each channel (R, G, B)
    retinexResult[0] = new Array();
    retinexResult[1] = new Array();
    retinexResult[2] = new Array();
    //for each line
    for (var i = 0; i < height; i++) {
        retinexResult[0][i] = new Array();
        retinexResult[1][i] = new Array();
        retinexResult[2][i] = new Array();
        for (var j = 0; j < width; j++) {
            retinexResult[0][i][j] = 0;
            retinexResult[1][i][j] = 0;
            retinexResult[2][i][j] = 0;
        }
    }

    //for all requested scales (MSR - multi scale retinex)
    for (var s = 0; s < sigmas.length; s++) {
        //gaussian blur coefs
        var coefs = GetGaussianCoeficients(sigmas[s])

        //Gaussian blur horizontal
        for (var i = 0; i < height; i++) {
            GaussianSmooth(imageDataData, bluredImageData, i * width * 4 + 0, 4, width, coefs);
            GaussianSmooth(imageDataData, bluredImageData, i * width * 4 + 1, 4, width, coefs);
            GaussianSmooth(imageDataData, bluredImageData, i * width * 4 + 2, 4, width, coefs);
        }

        //Gaussian blur vertical
        for (var i = 0; i < width; i++) {
            GaussianSmooth(bluredImageData, bluredImageData, i * 4 + 0, 4 * width, height, coefs);
            GaussianSmooth(bluredImageData, bluredImageData, i * 4 + 1, 4 * width, height, coefs);
            GaussianSmooth(bluredImageData, bluredImageData, i * 4 + 2, 4 * width, height, coefs);
        }

        //Retinex (SSR - single scale retinex)
        for (var i = 0; i < height; i++) {
            for (var j = 0; j < width; j++) {
                var index = (i * width + j) * 4;
                retinexResult[0][i][j] += weights[s] * Math.log((imageDataData[index + 0] + 1.0) / (bluredImageData[index + 0] + 1.0));
                retinexResult[1][i][j] += weights[s] * Math.log((imageDataData[index + 1] + 1.0) / (bluredImageData[index + 1] + 1.0));
                retinexResult[2][i][j] += weights[s] * Math.log((imageDataData[index + 2] + 1.0) / (bluredImageData[index + 2] + 1.0));
            }
        }

    }

    //Normalize result
    var min = retinexResult[0][0][0];
    var max = min;
    for (var i = 0; i < height; i++) {
        for (var j = 0; j < width; j++) {
            if (min > retinexResult[0][i][j]) min = retinexResult[0][i][j];
            else if (max < retinexResult[0][i][j]) max = retinexResult[0][i][j];
            if (min > retinexResult[1][i][j]) min = retinexResult[1][i][j];
            else if (max < retinexResult[1][i][j]) max = retinexResult[1][i][j];
            if (min > retinexResult[2][i][j]) min = retinexResult[2][i][j];
            else if (max < retinexResult[2][i][j]) max = retinexResult[2][i][j];
        }
    }

    //Put result to the new image
    var newImageData = ctx.createImageData(width, height);
    var newImageDataData = newImageData.data;
    var range = (max - min) / 255.0;
    for (var i = 0; i < height; i++) {
        for (var j = 0; j < width; j++) {
            //color restoration (MSRCR - multiscale retinex with color restoration)            
            var index = (i * width + j) * 4;
            var sum = imageDataData[index + 0] + imageDataData[index + 1] + imageDataData[index + 2] + 3.0;

            value = retinexResult[0][i][j] * beta * Math.log(alfa * (imageDataData[index + 0] + 1.0) / sum);
            value = (value - min) / range;
            if (value > 255) newImageDataData[index + 0] = 255;
            else if (value < 0) newImageDataData[index + 0] = 0;
            else
                newImageDataData[index + 0] = value;

            value = retinexResult[1][i][j] * beta * Math.log(alfa * (imageDataData[index + 1] + 1.0) / sum);
            value = (value - min) / range;
            if (value > 255) newImageDataData[index + 1] = 255;
            else if (value < 0) newImageDataData[index + 1] = 0;
            else
                newImageDataData[index + 1] = value;

            value = retinexResult[2][i][j] * beta * Math.log(alfa * (imageDataData[index + 2] + 1.0) / sum);
            value = (value - min) / range;
            if (value > 255) newImageDataData[index + 2] = 255;
            else if (value < 0) newImageDataData[index + 2] = 0;
            else
                newImageDataData[index + 2] = value;
        }

    }

    //set as no trasparent
    for (var i = 0; i < height; i++) {
        for (var j = 0; j < width; j++) {
            var index = (i * width + j) * 4;
            newImageDataData[index + 3] = 255;
        }
    }

    // copy the image data back onto the canvas
    ctx.putImageData(newImageData, 0, 0);
}

/* load image pointed by param_file to canvas */

function loadImage(imgSrc) {
    var canvas = document.getElementById("canvas");
    var ctx = canvas.getContext("2d");

    // add file:// if user specified local path
    if (imgSrc.indexOf("//") == -1 && imgSrc.indexOf(".") != 0) {
        imgSrc = "file:///" + imgSrc;
    }

    // load file into canvas
    var img = new Image();
    img.onload = function() {
        var width = img.width;
        var height = img.height;
        canvas.width = width;
        canvas.height = height;
        ctx.drawImage(img, 0, 0);
        // replace transparent with white
        try {
            imageData = ctx.getImageData(0, 0, width, height);
        } catch (e) {
            netscape.security.PrivilegeManager.enablePrivilege("UniversalBrowserRead");
            imageData = ctx.getImageData(0, 0, width, height);
        }
        for (var i = 0; i < height; i++) {
            for (var j = 0; j < width; j++) {
                index = (i * width + j) * 4;
                if (imageData.data[index + 3] == 0) {
                    imageData.data[index + 3] = 255;
                    imageData.data[index + 0] = 255;
                    imageData.data[index + 1] = 255;
                    imageData.data[index + 2] = 255;
                }
            }
        }
    }
    img.src = imgSrc;
}