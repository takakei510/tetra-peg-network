# 代表占有率の選定根拠

[Colabで開く](https://colab.research.google.com/github/takakei510/tetra-peg-network/blob/main/notebooks/Selection_basis_colab.ipynb)

「ランタイム → すべてのセルを実行」で実行します。Driveへのアクセスを許可すると、`MyDrive/tetra-peg-network/h2_report/selection_basis_colab` にPNG/SVG各2枚と選定結果・実行記録を保存します。データを手で配置する必要はありません。

入力は既存18,300試料・183条件の小さい集計CSV `summary.csv` です。各点は100試料平均とSEで、Pmaxと最大クラスタを1個除いたサイズ加重平均クラスタサイズSfinを描きます。プロットの元の平均・SEはローカルの18,300試料から再計算し一致を確認しています。新しいシミュレーションや原結合からの再集計ではありません。

選定規則は、ピークより低い側の平均Pmax≤5%の最大占有率、平均Sfin最大の測定点、ピークより高い側の平均Pmax≥80%の最小占有率です。L=128では0.41,0.43,0.60です。5%・80%は操作的な比較基準で相境界ではなく、ピークは厳密な転移点ではありません。既存データを見て決めた記述的な選定で、H2曲線の値を選定には使っていません。横軸は目標占有率です。

コードは `scripts/analysis/plot_selection_basis.py`。ノートブックのREVISIONでコード・CSVを同じコミットに固定しています。更新時はREVISIONも変更してください。実行時にSHA-256、版、条件を記録します。

このノートはH2_report_colab.ipynbとは別です。既存のH2ノートの描画条件は変更しません。
