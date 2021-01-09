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
2. (pipenv を使用する場合)pipenv の初期化を行う
   ```console
   $ pipenv install
   ```
3. seqtools を仮想環境にインストールする
   ```console
   $ pipenv run python seqtools/setup.py develop
   ```
4. `setting.yml`を生成する
   ```console
   $ pipenv run python src/initialize.py
   ```
   email は自分のものを記載する(Entrez.email に使う)


以降、自分の行いたい作業を行ってください

* csvを元にgbkを取得、`focus_rank`にしたがってディレクトリに分ける
```console
$ pipenv run python -m src.data_preparation.fetch
```
    * これをそれぞれ単品で行う場合は次の通り
    * csvを元にgbkを取得
```console
$ pipenv run get_gbk_by_csv
or 
$ pipenv run python -m src.data_preparation.fetch --source_csv <csvへのpath> --destination <gbkを保存するディレクトリ>
```
    * あるディレクトリに 
<!--  -->
<!-- - [organelles db](https://www.ncbi.nlm.nih.gov/genome/browse#!/organelles/)から取得した csv を元に genbank から gbk ファイルを取得する -->
<!--  -->
<!--   ```console -->
<!--   $ pipenv run python src/fetch_data/fetch.py [-i input_csv] [-d destination] -->
<!--   ``` -->
<!--  -->
<!--   `input_csv`, `destionation`の指定がないときは`setting.yml`に記載されたパスを使用する -->
<!--   `invalid_creatures`に記載されたパスに雑種などの使用しない生物を json 形式で書き出す -->
<!--  -->
<!-- - 取得した genbank ファイルの生物の分類学情報を取得する -->
<!--  -->
<!--   ```console -->
<!--   $ pipenv run python src/fetch_data/fetch_class.py [-i inputs...] [-d destination] -->
<!--   ``` -->
<!--  -->
<!--   `inputs`, `destionation`の指定がないときは`setting.yml`に記載されたパスを使用する -->
<!--  -->
<!-- - 取得した分類学情報をもとに gbk ファイルを`focus_rank`別にディレクトリわけする -->
<!--  -->
<!--   ```console -->
<!--   $ pipenv run python src/classification/classification.py [-i input_dir] -->
<!--   ``` -->
<!--  -->
<!-- - 取得した生物のゲノム配列を画像化する -->
<!--   ```console -->
<!--   $ pipenv run python src/generate_graph_img/generate_img.py [-w weight] [-wo weight_output] [-i inputfiles] [-d destination] -->
<!--   ``` -->
<!--   指定がないものは`setting.yml`から読み込む -->
