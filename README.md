# ミトコンドリアゲノム配列グラフ画像に基づく機械学習を用いた生物種の分類

## イントロダクション

一般に複数種の生物の進化的関係の評価方法としてアライメントが用いられているがこれは使用する生物種数$m$の遺伝配列がそれぞれ$n$の長さを持つとき、$\mathcal{O}{n^m}$の計算量を要する。

これは、タンパク質配列などの短い配列の比較であれば良いが、核ゲノム配列などの大規模な配列長の配列になると現実的な計算量とは言えない。
そこで、機械学習を用いて、短い時間でアライメントと同様の比較を行うというのが本研究の趣旨である。
また、機械学習に用いるデータとして、ミトコンドリアのゲノム配列を一定の規則により画像化したものを入力とし、入力された生物がどの綱に属するかを判別するのが現時点での目標である。

## 使い方
1. このリポジトリをクローンする。
    ```console
    $ git clone https://github.com/Omochice/master_thesis.git
    ```
2. (pipenvを使用する場合)pipenvの初期化を行う
    ```console
    $ pipenv sync
    ```
3. `setting.yml`を生成する
   ```console
   $ pipenv run python src/initialize.py
   ```
   emailは自分のものを記載する(Entrez.emailに使う)

* [organelles db](https://www.ncbi.nlm.nih.gov/genome/browse#!/organelles/)から取得したcsvを元にgenbankからgbkファイルを取得する
  ```console
  $ pipenv run python src/fetch_data/fetch.py [-i input_csv] [-d destination]
  ```
  `input_csv`, `destionation`の指定がないときは`setting.yml`に記載されたパスを使用する


* 取得したgenbankファイルの生物の分類学情報を取得する
  ```console
  $ pipenv run python src/fetch_data/fetch_class.py [-i inputs...] [-d destination]
  ```
  `inputs`, `destionation`の指定がないときは`setting.yml`に記載されたパスを使用する

* 取得した生物のゲノム配列を画像化する
  ```console
  $ pipenv run python src/generate_graph_img/generate_img.py [-w weight] [-wo weight_output] [-i inputfiles] [-d destination]
  ```
  指定がないものは`setting.yml`から読み込む