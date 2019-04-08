# 熱力学モデル

## 概要
このリポジトリにあるMAFit.cppは、熱力学モデルをC++で実装したものであり、研究用に作成したものです。  
両親媒性分子であるリン脂質は、単独で存在するモノマーの状態と、会合体として存在するミセルの状態の二つを取ります。
この熱力学モデルは、ある量のDHPC(1,2-dihexanoyl-sn-glycero-3-phosphocholine)を水に溶かしたときに、モノマーの量と、ミセルの量がどれだけかを見積もるのが目的です。
比熱を測る装置であるDSC(Differential Scanning Calorimetry)装置から得られた実験データをモデルにあてはめ、一番誤差関数が小さくなるようなパラメータを求めます。

## リポジトリの内容  
- MAFit.cpp 熱力学モデルをC++で実装したもの  
- MAFit     MAFit.cppをコンパイルしたバイナリーファイル(FreeBSD 10.3-STABLE FreeBSD clang version 3.4.1)  
- MAModel_short.pdf 熱力学モデルの詳細  
- 15mM, 16mM, 18mM, 20mM それぞれDHPC15mM, 16mM, 18mM, 20mMのDSCの実測データ(mMは濃度の単位で、mM=mmol/lです)  
- monomer, micelle それぞれモノマーとミセルの熱容量(J/K mol)のデータです  

## 使い方
コンパイルはclang++で  
<br/>
$ clang++ MAFit.cpp -o MAFit  
<br/>
とします。  
<br/>
$ MAFit -h  
<br/>
とすると、使い方を表示させることができます。  
主な使い方は次の通りです。  

- パラメーターを入力すると、それに対応したDSCの実測データを出力します。  
例えばDHPC15mMのDSC測定データを計算する場合は次のようにします。  
<br/>
$ MAFit -M 15 -n 20 -Tr 315 -b 30 -Cp1 monomer -Cp2 micelle  
<br/>  
計算の詳細を表示した後に、エンターを押すと、標準出力に計算結果を表示します。これをファイルに出力し、gnuplotなどのプロット用ソフトを用いてプロットします。  
<br/>  
$ MAFit -M 15 -n 20 -Tr 315 -b 30 -Cp1 monomer -Cp2 micelle > 15mM_model  
$ gnuplot  
$ plot "15mM_model" using 1:10  
<br/>  
- DSCの実測データと、パラメーターを動かす範囲を入力すると、パラメーターを動かしながら、誤差関数を計算して出力します。最後に誤差関数が一番小さくなるようなパラメーターを表示します。  

## モデル
以下のようなモデルを実装しています。
![image](./MAModel_short.png)


## DSCファイルのフォーマット
DSCの実験装置からの出力は次のようなフォーマットになっています。  

'# 温度(Celsius) 比熱(J/K)'  
0.000   0.00004792  
0.032   0.00002448  
0.069   0.00002077  
0.109   0.00002613  
  .         .  
  .         .  
  .         .  
  
