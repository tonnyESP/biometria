//////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                              //
//                             Simple image processing filters app                              //
//                               @author: Tonny Ruiz-Gijón                                      //
//                              @info: Biometria. MIARFID 2016                                  //
//                                                                                              //
//////////////////////////////////////////////////////////////////////////////////////////////////
// Caras de http://www.inf.unideb.hu/ipgd/FaceColorSegmentation/


// "Clase" ImageProcessing para el tratamiento de imagen en el navegaodr
function ImageProcessing(canvasName)
{

  /* ------------- Parámetros del canvas e imagen seleccionada -------------- */
  this.canvasName = canvasName;
  this.c   = document.getElementById(canvasName);
  this.ctx = this.c.getContext("2d");
  this.img = null;

  /* ---------- Canvas temporal donde realizar las convoluciones  ----------- */
  this.tmpc   = document.getElementById(canvasName+'_temporal');
  this.tmpctx = this.tmpc.getContext("2d");


  /* ------- Histograma y su función de distribución acumulada -------------- */
  this.histogram = Array(256).fill(0);
  this.histogCdf = Array(256).fill(0);

  /* ------------------- Media y varianza de la imagen ---------------------- */
  this.mean     = 0;
  this.variance = 0;

  /* --------------- Imagen integral que será computada --------------------- */
  this.integralImage = null;
  this.iic   = document.getElementById(canvasName+'_integralImage');
  this.iictx = this.iic.getContext("2d");



  // Se carga la imagen
  this.LoadImage();
};

// Crea el canvas temporal para las convoluciones de un width y height determinado
ImageProcessing.prototype.CreateTmpImageData = function(w, h)
{
  this.tmpctx.createImageData(w, h);

 // Se ajusta el tamaño del canvas al tamaño de la imagen
  $("#"+this.canvasName+"_temporal").attr("width", w);
  $("#"+this.canvasName+"_temporal").attr("height", h)
}

// Carga la imagen en en canvas y lo ajusta
ImageProcessing.prototype.LoadImage = function()
{
  // Se actualiza la imagen seleccionada
  this.img = document.getElementsByClassName("selected_image")[0];


  // Se ajusta el tamaño del canvas al tamaño de la imagen
  $("#"+this.canvasName).attr("width", this.img.width);
  $("#"+this.canvasName).attr("height", this.img.height)

  // Se ajusta el tamaño del canvas de la imgen integral al tamaño de la imagen
  $("#"+this.canvasName+'_integralImage').attr("width", this.img.width);
  $("#"+this.canvasName+'_integralImage').attr("height", this.img.height)

  // Se pinta la imagen seleccionada en el canvas
  this.ctx.drawImage(this.img, 0, 0, this.img.width, this.img.height);

  this.iictx.fillRect(0, 0, this.img.width, this.img.height);

  this.ComputeHistogram();
  this.ComputeMeanAndVariance();
  this.ComputeIntegralImage();
  // Se ajusta el tamaño de la gráfica con el histograma
  //$("#histogramChart").css("height", $(".dg").height());
}

// Devolvemos el array de pixels rgba
ImageProcessing.prototype.GetPixels = function() 
{
  return this.ctx.getImageData(0, 0, this.c.width, this.c.height);
};

// Devolvemos el array de pixels rgba del canvas temporal
ImageProcessing.prototype.GetTmpPixels = function() 
{
  return this.tmpctx.getImageData(0, 0, this.tmpc.width, this.tmpc.height);
};

// Devolvemos el array de pixels rgba del canvas de la imagen integral
ImageProcessing.prototype.GetIIPixels = function() 
{
  return this.iictx.getImageData(0, 0, this.iic.width, this.iic.height);
};
// Actualiza la información del canvas
ImageProcessing.prototype.PutPixels = function(imgData)
{
  // Actualizamos el canvas
  this.ctx.putImageData( imgData, 0, 0);

  // Computamos el histograma cada vez qeu apliquemos un nuevo filtro
  this.ComputeHistogram();
  // También se recalcula la media y varianza
  this.ComputeMeanAndVariance();
  // Se actualiza la imagen integral
  this.ComputeIntegralImage();

}
// Actualiza la información del canvas para la imagen integral
ImageProcessing.prototype.PutIIPixels = function(imgData)
{
  // Actualizamos el canvas
  this.iictx.putImageData( imgData, 0, 0);
}
// Devuelve el valor de rojo de un pixel representando la imagen en matrix 2D mediante un x,y
ImageProcessing.prototype.GetPixelValue = function(x, y)
{
  var imgData = this.GetPixels();

  var index = (x + y * imgData.width) * 4;
  var r = imgData.data[index+0];
  var g = imgData.data[index+1];
  var b = imgData.data[index+2];
  var a = imgData.data[index+3];

  return r;
}
// Asigna los valores de rgba a la imagen 1D
ImageProcessing.prototype.SetPixelValue = function(imgData, x, y, r, g, b, a) 
{
  g = b = r;
  a = 255;
  var index = (x + y * imgData.width) * 4;
  imgData.data[index+0] = r;
  imgData.data[index+1] = g;
  imgData.data[index+2] = b;
  imgData.data[index+3] = a;

}

// Ajusta a escala de grises aplicando:
// bnvalue = 0.299×r + 0.587×g + 0.114×b
ImageProcessing.prototype.BlackAndWhite = function(otherImageData) 
{
  if(otherImageData != undefined)
    var imgData = otherImageData;
  else
    var imgData = this.GetPixels();

  for (var i = 0; i < imgData.data.length; i += 4)
  {

    var bnvalue = Math.round( imgData.data[i] * 0.299 + imgData.data[i+1] * 0.587 + imgData.data[i+2] * 0.114);

    imgData.data[i] = imgData.data[i+1] = imgData.data[i+2] = bnvalue;
  }

  // Si es para devolver, damos la info
  if(otherImageData != undefined)
    return imgData;
  // Se aplica
  else
    this.PutPixels(imgData);
}

// Invierte los colores:
ImageProcessing.prototype.Invert = function() 
{
  var imgData = this.GetPixels();

  for (var i = 0; i < imgData.data.length; i += 4)
  {

    var bnvalue = Math.round( imgData.data[i] * 0.299 + imgData.data[i+1] * 0.587 + imgData.data[i+2] * 0.114);

    imgData.data[i]   = 255 - imgData.data[i];
    imgData.data[i+1] = 255 - imgData.data[i+1];
    imgData.data[i+2] = 255 - imgData.data[i+2];
    imgData.data[i+3] = 255;
  }

  // Se aplica
  this.PutPixels(imgData);
}

/* Al final se utiliza el kernel de 3x3 precompuesto
// Devuelve la gaussiana para un x, mu y sigma dado
function gaussian (x, mu, sigma) {
  return Math.exp( -(((x-mu)/(sigma))*((x-mu)/(sigma))) * 0.5 );
}
// Aplica gaussian Blur al canvas
ImageProcessing.prototype.GaussianBlur = function(kernelRadius)
{
  var sigma = this.mean;//kernelRadius * 0.5;
  var kernel = Array(2*kernelRadius+1).fill(Array(2*kernelRadius+1));
  var sum = 0;
  // compute values
  for (var row = 0; row < kernel.length; row++)
  {
    for (var col = 0; col < kernel[row].length; col++) {
      var x = gaussian(row, kernelRadius, sigma) * gaussian(col, kernelRadius, sigma);
      kernel[row][col] = x;
      sum += x;
    }
  }
  // Se debería normalizar
  //for (var row = 0; row < kernel.length; row++)
  //  for (var col = 0; col < kernel.length; col++)
  //    kernel[row][col] /= sum;

  // Se pasa de 2D a 1D
  var filter = [];
  for(var i = 0; i < kernel.length; i++)
  {
      filter = filter.concat(kernel[i]);
  }

  this.Convolute(filter, true);
}
*/

ImageProcessing.prototype.HistogramEqualization = function() 
{
    // Para ecualizar 
    var imgData = this.GetPixels();
    var norm = 255 / (imgData.width * imgData.height);

    for (var i = 0; i < imgData.data.length; i += 4) 
    {
        imgData.data[i] = this.histogram[imgData.data[i]] * norm;
    }

    this.PutPixels(imgData);
}

// Normalización global de la imagen
ImageProcessing.prototype.GlobalNormalization = function() 
{
  var imgData = this.GetPixels();

  this.ComputeMeanAndVariance();

  for (var i = 0; i < imgData.data.length; i += 4)
  {
    imgData.data[i]   = Math.round(imgData.data[i] - this.mean) / this.variance;
    imgData.data[i+1] = Math.round(imgData.data[i+1] - this.mean) / this.variance;
    imgData.data[i+2] = Math.round(imgData.data[i+2] - this.mean) / this.variance;
    imgData.data[i+3] = 255;   
  }

  this.PutPixels(imgData);
}

// Normalización Local
ImageProcessing.prototype.LocalNormalization = function()
{
  
  // Clon de la matriz de píxeles
  var pixels = this.GetPixels();
  // Se aplica el filtro gaussiano
  var means = this.GaussianFilter(sigma_mean, pixels);

  // Clonamos la matrix de píxeles
  var centeredPixels = this.GetPixels().slice(0);

  for (var i = 0; i < centeredPixels.data.length; i += 4)
    centeredPixels.data[i] = centeredPixels.data[i] - means[i];

  stds = np.sqrt(ndimage.gaussian_filter(im_centered*2, sigma_std))

}


// Convoluciona un filtro a la imagen. Se indica el filtro, si emplea el alfa, y si ha de devolver los valores o aplicarlos directamente
ImageProcessing.prototype.Convolute = function(weights, opaque, otherImageData) 
{
  // Tamaño del filtro y tamaño medio
  var side = Math.round(Math.sqrt(weights.length));
  var halfSide = Math.floor(side/2);

  // Info de la imagen
  var imgData   = this.GetPixels();
  var imgWidth  = imgData.width;
  var imgHeight = imgData.height;

  // Padding alrededor de la matriz de convolución
  var w = imgWidth;
  var h = imgHeight;

  // Se crea el canvas temporal si se ha de aplicar o Float32Array si se ha de devolver
  if(otherImageData != undefined)
    var toApply = { width: w, height: h, data: new Float32Array(w * h * 4) };    // RGBA
  else
    var toApply = this.tmpctx.createImageData(w, h);

  var dst = toApply.data;

  // Se recorren 
  var alphaFac = opaque ? 1 : 0;
  for (var y=0; y<h; y++) {
    for (var x=0; x<w; x++) {
      var sy = y;
      var sx = x;
      var dstOff = (y*w+x)*4;

      // Se calcula la suma ponderada de la imagen orignial dentro de la matriz de convolución
      var r = 0, 
          g = 0, 
          b = 0, 
          a = 0;

      for (var cy = 0; cy < side; cy++) 
      {
        for (var cx = 0; cx < side; cx++) 
        {
          var scy = sy + cy - halfSide;
          var scx = sx + cx - halfSide;
          if (scy >= 0 && scy < imgHeight && scx >= 0 && scx < imgWidth) 
          {
            var srcOff = (scy*imgWidth+scx)*4;
            var wt = weights[cy*side+cx];
            r += imgData.data[srcOff] * wt;
            g += imgData.data[srcOff+1] * wt;
            b += imgData.data[srcOff+2] * wt;
            a += imgData.data[srcOff+3] * wt;
          }
        }
      }
      dst[dstOff]   = r;
      dst[dstOff+1] = g;
      dst[dstOff+2] = b;
      dst[dstOff+3] = a + alphaFac * (255 - a);
    }
  }


  if(otherImageData != undefined)
    return toApply;
  // Se aplica
  else
    this.PutPixels(toApply);

};
// Aplica el kernel gaussiano de 3x3 
ImageProcessing.prototype.GaussianBlur = function(otherImageData)
{
  return this.Convolute( 
    [ 1/16, 1/8, 1/16, 
      1/8, 1/4, 1/8, 
      1/16, 1/8, 1/16,  ]
     , false, otherImageData)
}
// Aplica la convolucion para enfocar la imagen
ImageProcessing.prototype.Sharpen = function()
{
  this.Convolute( 
    [  0, -1,  0,
      -1,  5, -1,
       0, -1,  0 ]
     )
}

ImageProcessing.prototype.Sobel = function()
{
  var otherImageData = this.GetPixels();
  // se pasa a blanco y negro
  otherImageData = this.BlackAndWhite(otherImageData);

  var vertical = this.Convolute(
    [-1,-2,-1,
      0, 0, 0,
      1, 2, 1], true, otherImageData);
  var horizontal = this.Convolute(
      [-1,0,1,
     -2,0,2,
     -1,0,1], true, vertical);

  var toApply = this.tmpctx.createImageData(vertical.width, vertical.height);

  for (var i=0; i<toApply.data.length; i+=4) {
    var v = Math.abs(vertical.data[i]);
    toApply.data[i] = v;
    var h = Math.abs(horizontal.data[i]);
    toApply.data[i+1] = h
    toApply.data[i+2] = (v+h)/4;
    toApply.data[i+3] = 255;
  }
  
  this.PutPixels(toApply);
}
// Código original C https://github.com/HomeOfVapourSynthEvolution/VapourSynth-Retinex
//sigmaS=[25,80,250], [1,1,1] 0.001, 0.001
ImageProcessing.prototype.Retinex = function(sigmas, weights, alfa, beta)
{
  var width  = this.c.width;
  var height = this.c.height;

  var imgData = this.GetPixels();

  // Clonamos
  var bluredImgData  = this.ctx.createImageData(imgData.width, imgData.height);
  bluredImgData.data = imgData.data;
  // Se necesita crear la imagen con blur
  //bluredImgData = this.GaussianBlur(bluredImgData);


  // Lo que devolverá retinex, matriz 2D para cada canal
  var retinexResult = new Array();
    // Para R, G, B, alfa no tenido en cuenta
    retinexResult[0] = new Array();
    retinexResult[1] = new Array();
    retinexResult[2] = new Array();
  
  // Se inizaliza a 0
  for (var i = 0; i < height; i++) 
  {
      retinexResult[0][i] = new Array();
      retinexResult[1][i] = new Array();
      retinexResult[2][i] = new Array();

      for (var j = 0; j < width; j++) 
      {
          retinexResult[0][i][j] = 0;
          retinexResult[1][i][j] = 0;
          retinexResult[2][i][j] = 0;
      }
  }

  // Para todas las escalas (MSR - multi scale retinex)
  for (var s = 0; s < sigmas.length; s++) 
  {
      // Coeficientes para cada uno de los sigmas dados
      var coefs = GetGaussianCoeficients(sigmas[s])

      // Gaussian blur en horizontal
      for (var i = 0; i < height; i++) 
      {
          GaussianSmooth(imgData.data, bluredImgData, i * width * 4 + 0, 4, width, coefs);
          GaussianSmooth(imgData.data, bluredImgData, i * width * 4 + 1, 4, width, coefs);
          GaussianSmooth(imgData.data, bluredImgData, i * width * 4 + 2, 4, width, coefs);
      }

      // Gaussian blur en vertical
      for (var i = 0; i < width; i++) {
          GaussianSmooth(bluredImgData, bluredImgData, i * 4 + 0, 4 * width, height, coefs);
          GaussianSmooth(bluredImgData, bluredImgData, i * 4 + 1, 4 * width, height, coefs);
          GaussianSmooth(bluredImgData, bluredImgData, i * 4 + 2, 4 * width, height, coefs);
      }

      //Retinex (SSR - single scale retinex)
      for (var i = 0; i < height; i++) {
          for (var j = 0; j < width; j++) {
              var index = (i * width + j) * 4;
              retinexResult[0][i][j] += weights[s] * Math.log((imgData.data[index + 0] + 1.0) / (bluredImgData[index + 0] + 1.0));
              retinexResult[1][i][j] += weights[s] * Math.log((imgData.data[index + 1] + 1.0) / (bluredImgData[index + 1] + 1.0));
              retinexResult[2][i][j] += weights[s] * Math.log((imgData.data[index + 2] + 1.0) / (bluredImgData[index + 2] + 1.0));
          }
      }

  }

  // Se normaliza
  var min = retinexResult[0][0][0];
  var max = min;
  for (var i = 0; i < height; i++) 
  {
      for (var j = 0; j < width; j++) 
      {
          if (min > retinexResult[0][i][j]) 
            min = retinexResult[0][i][j];
          else if (max < retinexResult[0][i][j]) 
            max = retinexResult[0][i][j];

          if (min > retinexResult[1][i][j]) 
            min = retinexResult[1][i][j];
          else if (max < retinexResult[1][i][j]) 
            max = retinexResult[1][i][j];

          if (min > retinexResult[2][i][j]) 
            min = retinexResult[2][i][j];
          else if (max < retinexResult[2][i][j]) 
            max = retinexResult[2][i][j];
      }
  }

  // Guardamos el resultado en una nueva imagen
  var newImageData = this.ctx.createImageData(width, height);
  var newImageDataData = newImageData.data;
  var range = (max - min) / 255.0;
  for (var i = 0; i < height; i++) {
      for (var j = 0; j < width; j++) {
          //color restoration (MSRCR - multiscale retinex with color restoration)            
          var index = (i * width + j) * 4;
          var sum = imgData.data[index + 0] + imgData.data[index + 1] + imgData.data[index + 2] + 3.0;

          value = retinexResult[0][i][j] * beta * Math.log(alfa * (imgData.data[index + 0] + 1.0) / sum);
          value = (value - min) / range;
          if (value > 255) newImageDataData[index + 0] = 255;
          else if (value < 0) newImageDataData[index + 0] = 0;
          else
              newImageDataData[index + 0] = value;

          value = retinexResult[1][i][j] * beta * Math.log(alfa * (imgData.data[index + 1] + 1.0) / sum);
          value = (value - min) / range;
          if (value > 255) newImageDataData[index + 1] = 255;
          else if (value < 0) newImageDataData[index + 1] = 0;
          else
              newImageDataData[index + 1] = value;

          value = retinexResult[2][i][j] * beta * Math.log(alfa * (imgData.data[index + 2] + 1.0) / sum);
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
  this.PutPixels(newImageData);
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

// Computa el histograma y la función de distribución acumulada de este
ImageProcessing.prototype.ComputeHistogram = function()
{

  var imgData = this.GetPixels();

  // Solo gastaremos el color rojo, porque en principio debe estar la imagen en blanco y negro
  for (var i = 0; i < imgData.data.length; i += 4)
    this.histogram[imgData.data[i]] ++;

  // Píxeles en la imagen (no es necesario que esté dividido por 4, que son los 4 valores de rgba. width y height devuelven el tamaño bueno)
  var n = (imgData.width * imgData.height) ;

  // Normalización histograma (flat)
  for (var i = 0; i < this.histogram.length; i++)
    this.histogram[i] = this.histogram[i] / n ;

  // Función de distribución acumulada
  this.histogCdf[0] = this.histogram[0];
  for (var i = 1; i < this.histogram.length; i++)
    this.histogCdf[i] = this.histogram[i] + this.histogCdf[i-1];

  // Ploteamos para ver qué aspecto tiene el histograma y su cdf
  PlotHistogram(this.histogram, this.histogCdf);
}

// Computa la media y la varianza en la imagen
ImageProcessing.prototype.ComputeMeanAndVariance = function()
{
  var data = this.GetPixels().data;

  var total = 0;
  var count = 0;
  for (var i = 0; i < data.length; i += 4){
    total += parseFloat(data[i]);
    count ++;
  };
  var mean = total / count;

  // reset counter
  count = 0;
  var totalVariance = 0;
  for (var i = 0; i < data.length; i += 4){
    var vari = data[i] - mean;
    totalVariance += (vari * vari);

    count++;
  };

  var variance = Math.sqrt( totalVariance / count);

  this.mean     = mean;
  this.variance = variance;
}

ImageProcessing.prototype.ComputeIntegralImage = function()
{
  var imgData = this.GetPixels();
  var pixels = imgData.data;
  var w = imgData.width;
  var h = imgData.height;

  var iiData = this.GetIIPixels();


  var imgIntegral = Array(w).fill(Array(h));

  // Integral en la posición actual
  var integral = 0;

  // Últimos valores de la imagen integral
  var lx, ly, lxy;

  for(var i = 0; i < w; i++)
  {
      for(var j = 0; j < h; j++) 
      {
        // índice en el array 1D de rgba
        var index = (i + j * imgData.width) * 4;
        // Suma de los valores de los píxeles
        // Comprobacióin borde inferior
        lx  = i > 0 ? imgIntegral[i-1][j] : 0;
        ly  = j > 0 ? imgIntegral[i][j-1] : 0;
        lxy = i > 0 && j>0 ? imgIntegral[i-1][j-1] : 0;

        imgIntegral[i][j] = pixels[index] + lx + ly - lxy; 

        var bnvalue = Math.round( imgIntegral[i][j] * 0.299 + imgIntegral[i][j] * 0.587 + imgIntegral[i][j] * 0.114);

        iiData.data[index+0] = bnvalue / 255;
        iiData.data[index+1] = bnvalue / 255;
        iiData.data[index+2] = bnvalue / 255;
        iiData.data[index+3] = 255;

      }
  }

  this.PutIIPixels(iiData);

}



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