# 熱力学モデル

## 概要
このリポジトリにあるMAFit.cppは、熱力学モデルをC++で実装したものであり、研究用に作成したものです。  
両親媒性分子であるリン脂質は、水中で単独で存在するモノマーの状態と会合体として存在するミセルの状態の二つを取ります。
この熱力学モデルは、ある量のリン脂質を水に溶かしたときに、モノマーの量と、ミセルの量がどれだけかを見積もるのが目的です。
比熱を測る装置であるDSC(Differential Scanning Calorimetry)装置から得られた実験データをモデルにあてはめ、一番誤差関数が小さくなるようなパラメータを求めます。


## 使い方
コンパイルはclang++で  
$ clang++ MAFit.cpp -o MAFit  
とします。-hをオプションにつけると、使い方を表示させることができます。  
主な使い方は次の通りです。  

- パラメーターを入力すると、それに対応したDSCの実測データを出力します。
- DSCの実測データと、パラメーターを動かす範囲を入力すると、パラメーターを動かしながら、誤差関数を計算して出力します。最後に誤差関数が一番小さくなるようなパラメーターを表示します。  

## 熱力学モデル
ΔH=ΔCp(T-Tr)

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
  
