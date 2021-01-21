# ミトコンドリアゲノム配列グラフ画像に基づく機械学習を用いた生物種の分類

## イントロダクション

一般に複数種の生物の進化的関係の評価方法としてアライメントが用いられているがこれは使用する生物種数$m$の遺伝配列がそれぞれ$n$の長さを持つとき、$\mathcal{O}{n^m}$の計算量を要する。

これは、タンパク質配列などの短い配列の比較であれば良いが、核ゲノム配列などの大規模な配列長の配列になると現実的な計算量とは言えない。
そこで、機械学習を用いて、短い時間でアライメントと同様の比較を行うというのが本研究の趣旨である。
また、機械学習に用いるデータとして、ミトコンドリアのゲノム配列を一定の規則により画像化したものを入力とし、入力された生物がどの綱に属するかを判別するのが現時点での目標である。

## 使い方

1. リポジトリをクローンする
```console
$ git clone https://github.com/Omochice/machine-genome-classification.git
```
2. pipenvを使って仮想環境を再現する
python3.8系を指定しているので別途pyenvなどで指定する
```console
$ cd machine-genome-classification
$ pipenv install 
```
3. 初期設定をする
```console
$ pipenv run python src/initialize.py
```
`setting.yml`が生成されるので`<FILL IN>`要素を環境に合わせて設定する

4. データのソースを[Genome List - Genome - NCBI](https://www.ncbi.nlm.nih.gov/genome/browse#!/organelles/)から取得する
Replicons列を見てデータの取得等をするので最低限その列は必要

5. gbkファイルを取得する
```console
$ pipenv run fetch
```

6. 画像を生成する
```console
$ pipenv run gen_img
```
