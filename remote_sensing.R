## --------------------------------------------------------------------
# title: "Rを用いたドローンマルチスペクトル画像解析"
# author: "株式会社エコリス 水谷貴行"
# date: "2018年7月30日"
## --------------------------------------------------------------------

#install.packages("tidyverse","sf","rgeos","raster","dismo",dependencies = T)
library(tidyverse)
library(raster)
library(sf)
library(ggplot2)
library(dismo)
library(rgeos)

## --------------------------------------------------------------------
# RedEdgeスペクトルカメラ画像処理
## --------------------------------------------------------------------

## --------------------------------------------------------------------
## 画像の読み込み、下処理、書き出し

#データを読み込みバンドを名を設定
rededge <- stack('flight20180703_rededge.tif')
rededge <- subset(rededge,1:5) #アルファバンドを除外
names(rededge) <- c('blue','green','red','NIR','rededge')
rededge

#切り取り(最初は手動で実施して値を取得)
#rededgecrop<-raster::select(rededge) #手動で選択する場合
e<-extent(24905,24992,2711,2824)
rededgecrop<-crop(rededge,e)

#書き出し
crs(rededgecrop)<-CRS("+init=epsg:2451")
writeRaster(rededgecrop,filename = "cropped-rededge.tif", format = "GTiff", overwrite = TRUE)

## --------------------------------------------------------------------
## 画像の読み込み、情報確認、表示

#データを読み込み、バンド名を指定、情報確認
rededge <- stack('cropped-rededge.tif')
names(rededge) <- c('blue','green','red','NIR','rededge')
crs(rededge)
ncell(rededge)
dim(rededge)
res(rededge)
nlayers(rededge)

#表示（blueバンド、TrueColor、FalseColor、NaturalColor）
plot(rededge,1,col=gray.colors(30))
plotRGB(rededge, r = 3, g = 2, b = 1, axes = TRUE, stretch = "lin", main = "RedEdge True Color Composite")
plotRGB(rededge, r = 4, g = 3, b = 2, axes = TRUE, stretch = "lin", main = "RedEdge False Color Composite")
plotRGB(rededge, r = 3, g = 4, b = 2, axes = TRUE, stretch = "lin", main = "RedEdge Natural Color Composite")

## --------------------------------------------------------------------
## バンド間の相関

#各バンドを抽出し、16bitの値から0～1の反射率に変換
b_red<-subset(rededge,3)/65536
b_nir<-subset(rededge,4)/65536
b_rededge<-rededge[[5]]/65536
red_nir<-stack(b_red,b_nir)
red_rededge<-stack(b_red,b_rededge)

#相関グラフ
pairs(red_nir, main = "Scatterplot between Red and NIR")
pairs(red_rededge, main = "Scatterplot between Red and Rededge")


## --------------------------------------------------------------------
## mapの使い方について



## --------------------------------------------------------------------
## スペクトルプロファイルの作成

#土地被覆の読み込み
landcover<-read_sf("landcover.shp",options=c("ENCODING=CP932"))
plot(landcover)
plot(landcover["class"],key.width = lcm(3))
plot(landcover["subclass"],key.width = lcm(3))

#サブクラスの一覧をベクトルで取得
subcl<-landcover %>% 
  st_set_geometry(NULL) %>% #sfからデータフレームに変換
  dplyr::distinct(subclass,.keep_all=FALSE) %>% #重複を削除 
  .$subclass #ベクトルを返す

#サブクラスごとにランダム点を300ポイント作成
subclass_sampling<-function(x){
  filter(landcover,subclass==x) %>% 
    st_sample(300)
}
pntlist<-subcl %>% map(subclass_sampling) 
pnt<-pntlist %>% reduce(c)

#ポイントの属性にclassとsubclassを追加
pnt<-st_intersection(landcover,pnt)
plot(pnt,pch=".",cex=3)

#各ポイントでのバンド値を抽出
band_df <- raster::extract(rededgecrop, as(pnt,"Spatial"),df=TRUE)
#ポイントの属性にバンド値（反射率に変換）を追加
pnt<-pnt %>% bind_cols(band_df/65536)

#土地被覆(大区分)ごとのバンドの平均値、最小値、最大値
df_class_mean<-pnt %>% 
  st_set_geometry(NULL) %>% 
  dplyr::select(-ID,-subclass) %>% 
  group_by(class) %>% 
  summarise_all(mean)%>% 
  tidyr::gather(key=Bands,value=Reflectance,-class)


df_class_min<-bind_cols(ptsamp,df) %>% st_set_geometry(NULL) %>% 
  dplyr::select(class,blue,green,red,rededge,NIR) %>% 
  group_by(class) %>% summarise_all(min)%>% 
  tidyr::gather(key=Bands,value=Ref_min,-class)

df_class_max<-bind_cols(ptsamp,df) %>% st_set_geometry(NULL) %>% 
  dplyr::select(class,blue,green,red,rededge,NIR) %>% 
  group_by(class) %>% summarise_all(max)%>% 
  tidyr::gather(key=Bands,value=Ref_max,-class)

df_class<-df_class_mean %>% left_join(df_class_min) %>% left_join(df_class_max) 

#土地被覆(小区分)ごとのバンドの平均値、最小値、最大値
df_subclass_mean<-bind_cols(ptsamp,df) %>% st_set_geometry(NULL) %>% 
  dplyr::select(subclass,blue,green,red,rededge,NIR) %>% 
  group_by(subclass) %>% summarise_all(mean)%>% 
  tidyr::gather(key=Bands,value=Reflectance,-subclass)

df_subclass_min<-bind_cols(ptsamp,df) %>% st_set_geometry(NULL) %>% 
  dplyr::select(subclass,blue,green,red,rededge,NIR) %>% 
  group_by(subclass) %>% summarise_all(min)%>% 
  tidyr::gather(key=Bands,value=Ref_min,-subclass)

df_subclass_max<-bind_cols(ptsamp,df) %>% st_set_geometry(NULL) %>% 
  dplyr::select(subclass,blue,green,red,rededge,NIR) %>% 
  group_by(subclass) %>% summarise_all(max)%>% 
  tidyr::gather(key=Bands,value=Ref_max,-subclass)

df_subclass<-df_subclass_mean %>% left_join(df_subclass_min) %>% left_join(df_subclass_max) 
df_subclass

df_class2 <- transform(df_class, Bands= factor(Bands, levels = c("blue", "green", "red", "rededge", "NIR")))

#大区分のバンドごとの反射率をグラフプロット
df_class2 %>% #dplyr::filter(class=="構造物") %>% 
ggplot(aes(x=Bands ,y=Reflectance,col=class))+
  #geom_ribbon(aes(ymin = Ref_min, ymax = Ref_max, group = class,fill=class), alpha = .2 )+
  geom_line(aes(group = class),size=1)+
  geom_point()+
#  coord_cartesian(ylim=c(0,0.4))+
  labs(title = "Spectral Profile from RedEdge")

df_subclass2 <- transform(df_subclass, Bands= factor(Bands, levels = c("blue", "green", "red", "rededge", "NIR")))

#小区分のバンドごとの反射率をグラフプロット
df_subclass2 %>% 
  ggplot(aes(x=Bands ,y=Reflectance,col=subclass))+
#  geom_ribbon(aes(ymin = Ref_min, ymax = Ref_max, group = subclass,fill=subclass), alpha = .2 )+
  geom_line(aes(group = subclass),size=1)+
  geom_point()+
#  coord_cartesian(ylim=c(0,0.4))+
  labs(title = "Spectral Profile from RedEdge")

#NDVIの計算
VI <- function(img, k, i) {
  bk <- img[[k]]
  bi <- img[[i]]
  vi <- (bk - bi) / (bk + bi)
  return(vi)
}

par(mfrow=c(1,2))
# For RedEdge NIR = 4, red = 3.
ndvi <- VI(rededgecrop,4, 3)
plot(ndvi, col = rev(terrain.colors(10)), main = 'RedEdge-NDVI')

# For RedEdge NDRE = 4, rededge = 5.
ndre <- VI(rededgecrop,4, 5)
plot(ndre, col = rev(terrain.colors(10)), main = 'RedEdge-NDRE')

# view histogram of data
hist(ndvi,
     main = "Distribution of NDVI values",
     xlab = "NDVI",
     ylab="Frequency",
     col = "wheat",
     xlim = c(-0.5, 1),
     breaks = 30,
     xaxt = 'n')
axis(side=1, at = seq(-0.5,1, 0.05), labels = seq(-0.5,1, 0.05))

veg <- calc(ndvi, function(x){x[x < 0.4] <- NA; return(x)})
plot(veg, main = 'Veg cover')

land <- reclassify(ndvi, c(-Inf,0.25,NA,0.25,0.3,1,0.3,Inf,NA))
plot(land, main = 'What is it?')
plotRGB(rededgecrop, r = 1, g = 2, b = 3, axes = TRUE, stretch = "lin", main = "Landsat False Color Composite")
plot(land, add = TRUE, legend = FALSE)

vegc <- reclassify(veg, c(-Inf,0.25,1, 0.25,0.3,2, 0.3,0.4,3, 0.4,0.5,4, 0.5,Inf, 5))
plot(vegc,col = rev(terrain.colors(4)), main = 'NDVI based thresholding')


#主成分分析
set.seed(1)
sr <- sampleRandom(rededgecrop, 10000)
plot(sr[,c(3,4)], main = "NIR-Red plot")
pca <- prcomp(sr, scale = TRUE)
pca
screeplot(pca)
pci <- predict(rededgecrop, pca, index = 1:2)
plot(pci[[1]])
pc2 <- reclassify(pci[[2]], c(-Inf,0,1,0,Inf,NA))
par(mfrow = c(1,2))
plotRGB(rededgecrop, r = 3, g = 2, b = 1, axes = TRUE, stretch = "lin", main = "Landsat False Color Composite")
plot(pc2, legend = FALSE, add = TRUE)


#教師なし分類(5バンド使用)
set.seed(99)
km <- kmeans(values(rededgecrop), centers = 6, iter.max = 500, nstart = 3, algorithm="Lloyd")
knr <-raster(rededgecrop,1)
knr[] <- km$cluster
par(mfrow = c(1,2))
plotRGB(rededgecrop, r = 3, g = 2, b = 1, axes = TRUE, stretch = "lin", main = "RedEdge True Color Composite")
plot(knr, main = 'Unsupervised classification',breaks=0:6,col=c("gray","orange","white","darkgreen","lightgreen","brown"))


#教師なし分類(ndvi使用)
ndvi[is.na(ndvi)]<-0
set.seed(99)
km <- kmeans(values(ndvi), centers = 5, iter.max = 500, nstart = 3, algorithm="Lloyd")
knr <- ndvi
knr[] <- km$cluster
par(mfrow = c(1,2))
plotRGB(rededgecrop, r = 3, g = 2, b = 1, axes = TRUE, stretch = "lin", main = "RedEdge True Color Composite")
plot(knr, main = 'Unsupervised classification',breaks=0:5,col=c("gray","green","brown","darkgreen","lightgreen"))

#教師あり分類
#カテゴリごとにランダム点を作成
cl<-landcover %>% st_set_geometry(NULL) %>% dplyr::distinct(subclass,.keep_all=FALSE) %>% .$subclass
pnt<-cl %>% map(~filter(landcover,subclass==.) %>% st_sample(300)) %>% reduce(c)
#土地被覆を属性に追加
ptsamp<-st_intersection(landcover,pnt)
plot(ptsamp["subclass"])

#バンドの値をランダム点で抽出
df <- raster::extract(rededgecrop, as(ptsamp,"Spatial"),df=TRUE)
#反射率には変換しない
#df <- df/65536
df <- df %>% dplyr::select(-ID)
#ランダム点とバンドの値を結合
sampdata<-bind_cols(ptsamp,df)
sampdata2<-sampdata %>% st_set_geometry(NULL) %>% dplyr::select(-class)
#トレーニングデータの作成
j <- kfold(sampdata2, k = 5, by = sampdata2$class)

training2011 <- sampdata2[j!= 1, ] # selected the rows where j equals 1
validation2011 <- sampdata2[j == 1, ] # selected the rows where j equals [2:k]

#training2011 <- sampdata2

#決定木で学習
library('rpart')
# Train the model
cart <- rpart(as.factor(subclass)~., data = training2011, method = 'class', minsplit = 5)
plot(cart, uniform=TRUE, main="Classification Tree")
text(cart, cex = 0.8)
#学習結果でラスタから予測
pr2011 <- predict(rededgecrop, cart, type='class', progress = 'text')
#予測結果をプロット
library(rasterVis)
classcolor <- c("gray", "black", "white", "darkgreen","brown","lightgreen","lightblue", "blue","orange","green","yellow")
#classcolor <- rainbow(11)
levelplot(pr2011, maxpixels = 1e6,
          col.regions = classcolor,
          scales=list(draw=FALSE),
          main = "Decision Tree classification of RedEdge camera")

###精度の評価
prclass <- predict(cart, validation2011[,2:ncol(validation2011)], type='class')
conmat <- data.frame(prediction = prclass, reference = validation2011$subclass)
conmat <- table(conmat)
print(conmat)
n <- sum(conmat)
print(n)
nc <- nrow(conmat)
print(nc)
diag <- diag(conmat)
rowsums <- apply(conmat, 1, sum)
colsums <- apply(conmat, 2, sum)
p <- rowsums / n
q <- colsums / n
OA <- sum(diag) / n
expAccuracy <- sum(p*q)
kappa <- (OA - expAccuracy) / (1 - expAccuracy)
print(OA)
print(kappa)
PA <- diag / colsums
UA <- diag / rowsums
outAcc <- data.frame(producerAccuracy = PA, userAccuracy = UA)
print(outAcc)

########
#SEQUOIA
########
sequoia <- stack('flight20180702_ortho.tif')
sequoia <- subset(sequoia,1:4)
names(sequoia) <- c('green','red','NIR','rededge')

#情報確認
crs(sequoia)
ncell(sequoia)
dim(sequoia)
res(sequoia)
nlayers(sequoia)
crs(sequoia)<-CRS('+init=EPSG:2451')

#表示
plot(sequoia,3,col=gray.colors(30))
tanbo<-read_sf("tanbo.shp",options=c("ENCODING=CP932"))
croped<-crop(sequoia,tanbo)
sequoia<-mask(croped,tanbo)
plot(sequoia)

#表示（切り取り範囲）
par(mfrow=c(1,2))
plotRGB(sequoia, r = 3, g = 2, b = 1, axes = TRUE, stretch = "lin", main = "RedEdge True Color Composite")
plotRGB(sequoia, r = 4, g = 3, b = 2, axes = TRUE, stretch = "lin", main = "RedEdge False Color Composite")
plotRGB(sequoia, r = 3, g = 4, b = 2, axes = TRUE, stretch = "lin", main = "RedEdge Natural Color Composite")

#書き出し
writeRaster(sequoia,filename = "cropped-sequoia.tif", format = "GTiff", overwrite = TRUE)

#ノイズ除去

nir<-raster(sequoia,3)
nir<- focal(nir, w=matrix(1/(3*3),nrow=3,ncol=3))
#ext<-select(nir)
ext<-extent(24790, 24792, 2761, 2763)
plot(nir,col=gray.colors(30),ext=ext)

gf <- focalWeight(nir, 0.1, "circle")
gf <- ifelse(gf == 0, 0, 1)
rg <- focal(nir, w=gf,fun=max)
plot(rg,ext=ext)
#rg<-crop(rg,ext)
pol <- rasterToPolygons(rg,n=4,dissolve=TRUE,fun=function(x){x>11000})
sfpol<-st_as_sf(pol) %>% st_cast("POLYGON")
plot(sfpol)
sfpol<-sfpol %>% mutate(AREA=st_area(.))
sfpol
fpol<-sfpol %>% filter(AREA>units::as_units(80,"cm^2"))
plot(fpol %>% st_geometry())
cent<-st_centroid(fpol)

plot(nir,col=gray.colors(30))
plot(fpol,add=T,col=NA)
plot(cent,add=T,pch=".",cex=2)

?st_make_grid
st_make_grid(n=rev(dimyx))

rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
p<-c(24787.8,2755.4)
g = st_make_grid(n=c(38,12), offset = p, cellsize = c(0.5,0.5))#
cp<-st_sfc(st_point(p))


g<-(g-cp)*rot(-pi*50/180)+cp
plot(nir,col=gray.colors(30))
plot(g,axes=T,add=T)
plot(cent,add=T,pch=".",cex=2)

#メッシュで集計
g<-st_set_crs(g,2451)
cnt<-g %>% st_intersects(cent,sparse=T) %>% map_int(length)
st_sf(g) %>% mutate(ct=cnt) %>% plot


#writeRaster(nir,filename = "sequoia-nir.tif", format = "GTiff", overwrite = TRUE)
#st_write(cent, "cent.shp", delete_layer = TRUE)

plot(nir,ext=ext,col=gray.colors(30))
plot(nir,ext=ext,col=gray.colors(30))

#NDVIの計算
VI <- function(img, k, i) {
  bk <- img[[k]]
  bi <- img[[i]]
  vi <- (bk - bi) / (bk + bi)
  return(vi)
}

#par(mfrow=c(1,2))

#NDVIの平均メッシュ表示
# For SEQUOIA NIR = 3, red = 2.
ndvi <- VI(sequoia,3, 2)
plot(ndvi, col = rev(terrain.colors(10)), main = 'SEQUOIA-NDVI')
plot(g,axes=T,add=T)
df<-raster::extract(ndvi,as(g,"Spatial"),fun=mean,df=TRUE)
st_sf(g) %>% mutate(ndvi=df$layer) %>% filter(ndvi>0.55) %>% plot(pal = rev(heat.colors(10)))

# For RedEdge NIR = 3, rededge = 4.
ndre <- VI(sequoia,3, 4)
plot(ndre, col = rev(terrain.colors(10)), main = 'SEQUOIA-NDRE')
