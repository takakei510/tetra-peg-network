# 基本クラスタ量とH₂解析のColabノート

## 基本クラスタ量（H₂を使わない）

[基本クラスタ量をColabで開く](https://colab.research.google.com/github/takakei510/tetra-peg-network/blob/main/notebooks/Cluster_basics_colab.ipynb)

「ランタイム → すべてのセルを実行」で、PC版と同じ7種類の図を作ります。入力は今回の183条件・各100試料の `samples.csv` です。初回は入力がDriveにない場合だけ、PCの結果フォルダからこのCSVをアップロードしてください。以降はDriveの `tetra-peg-network/cluster_basics/input/samples.csv` を再利用します。`summary.csv` や `size_relation.csv` は入力にできません。

出力はDriveの `tetra-peg-network/cluster_basics/figures_<日時>` です。PNG/SVG各7枚、集計CSV、定義・目的・説明例を記した `figure_guide.md`、入力とコードのSHA256を記した `run_metadata.json` が保存されます。最後のボタンで `cluster_basics_share.zip` をダウンロードし、ChatGPTに添付できます。コードと設定をコミットに固定し、入力は試料ごとの比を計算してから等重みで平均します。シミュレーションは実行しません。

# H₂ report notebook

[Open in Colab](https://colab.research.google.com/github/takakei510/tetra-peg-network/blob/main/notebooks/H2_report_colab.ipynb)

## Colabでの実行

1. 上のリンクからGitHubの最新版ノートブックを開きます。既に開いたDriveのコピーにはGitHubの変更は自動反映されません。
2. 上からセルを実行し、Google Driveへのアクセスを許可します。
3. データ設定セルのROOTを確認します。既定は `/content/drive/MyDrive/tetra-peg-network/h2_report` です。このフォルダには `size_relation.csv` が必要です。
4. 必要に応じてLとPHISを変更します。
5. 図6ファイル（PNG/SVG各3種類）と `run_metadata.json` がDriveの `figures_colab` に保存されます。同名出力は更新されます。

## コードの管理と再現

描画処理は `scripts/analysis/plot_report.py` で管理しています。Drive内の `plot_report.py` は不要になりましたが、削除する必要はありません。
ノートブックのCODE_REFは描画コードのコミット番号に固定しています。描画コードを変更した場合は、CODE_REFも更新してノートブックを保存してください。
`main` を指定する場合も、実行時に解決したコミット番号を表示し、記録します。取得に失敗した場合に古いDriveのコードへ切り替えることはありません。

実行記録にはコミット番号、描画条件、入力CSV・コードのSHA-256、Python・matplotlib版を保存します。入力CSV・原データ・出力画像はこの変更ではGitHubに追加していません。新たに実行した結果はColabからGitHubへ保存しない限りリポジトリには反映されません。

## 集計の意味

既存の描画処理をそのまま追加しています。100試料の既存集計CSV向けの注記、H₂割合の縦軸0〜11%、寄与試料数50を目安にする表示は元の仕様です。異なる試料数や範囲に使う場合は注記・軸も見直してください。
横軸はサイズ区間を等間隔に並べたカテゴリ軸です。クラスタ数は区間幅で割っていません。H₂割合は試料内でクラスタを等重み平均し、その区間のある試料間で等重み平均しています。

GitHub上のノートブックは古い図と新しいコードの混同を避けるため、今回の更新時に保存済み出力を消しています。実行すると図が再表示されます。
