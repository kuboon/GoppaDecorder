# GoppaDecorder(WhiteBlood)

１．プログラム実行に当たって(oplib.c)　20200304

make

とやれば出来ます。
出力の第一段目はエラーの位置と値をユークリッドアルゴリズムで計算したもの、(decode関数)
第二は復号法であるパターソンアルゴリズムでバイナリエラーの位置を出力したものになります。(pattarson関数)

20200313

とりあえず関数encryptで512ビット秘密鍵は暗号化復号化できるようになったので、Niederreiter暗号は完成していると言っていい。
残るのは鍵をバイナリにして生成と復号をするかというUIの部分だけだが、ここでPGPのような暗号アプリを作るつもりはないので、
今後コードを利用して必要なアプリや各自応用をしてほしい。
自分もできる限りコードの可読性を上げるつもりなのでメンテナンスは続けていく。
バイナリ符号を扱うことができるようにするところまでは最低限やっておきたい。

20200315

今windows版を作っています。こうして違う環境とデバッグをすると、一見うまく動いているように見えるだけという怪しいコードが結構ありますね。
できるだけエラーを減らします。それにしてもｇｃｃが性能がいいのかほとんどエラーは出ません。24時間耐久テストを実行してます。

20200320

AMD性最適化コンパイラ（https://developer.amd.com/amd-aocc/）

AMDを使っているならこちらで試してみてください。

このプログラムは今まで学習してきたことの集大成として実験的に古典的Goppa符号を複合するアルゴリズムを実装したものです。
応用としては公開鍵暗号やブロックチェーンがあります。
あくまでも自分用なのでおすすめはしません。
特徴としては、一変数多項式を単項式の配列で表現するという直感的な四則演算を設計に取り入れました。
アルゴリズムは設計者の思考パターンを反映するのかもしれません。
設計方法が悪いので最適化には向いてないかもしれませんが、とにかくコードの美しさよりまともに動くプログラムを優先的に作ろうとしたものです。
外部ライブラリを必要とせず、このリポジトリだけで完結しています。
とりあえず５１２ビットの秘密鍵を暗号化復号化ができるようになったので、今後はバイナリ符号に対応するのとメンテナンスとして可読性を上げていくようにします。
ちょっとづつですが毎日いらない変数を消したりしてます。
分割コンパイルするとオブジェクト間の最適が出来ないらしいので、意図的に分割コンパイルはしてません。
メインのソースにインクルードする形で一つのファイルにまとめています。
どちらかというとアプリではなく、特殊な科学計算なので速度優先で書いてます。

https://classic.mceliece.org/

こちらにプロが作った比較用のサンプルコードがあります。


２．作るに至った経緯

耐量子計算機暗号のコアの一つとして、古典的Goppa符号の復号法ができれば良い。ここでは暗号化に必要な関数を作りました。
現在のところ、2の拡大体に関する復号です。バグがあります。ベータバージョンです。マニュアルもまだ書いてません。
QC-MDPCなんか邪道だと思っているので、オーソドックスな実装にしました。実用に耐えるようなパラメータの実装はありませんでした。
あるにはあるのですが、外部ライブラリに頼った実装ばかりで、自己完結的に作成できるようなものがないので作りました。
現時点（20200306）では、バイナリ公開鍵には対応していません。

このコードを実装するのに使った書籍は、Oliver PretzelのError-Correcting Codes and Finite Fieldsです。

その他英語の論文についてはdiary.txtを参照してください。

このコードの最終目的は軍事レベルの暗号を作ることで、Niederreiter暗号を目指しています。鍵はでかいです。
毎月たくさんの新型が生まれている中で、なぜ改良版を使わないのかという人もいるかもしれませんが、
一番古いオーソドックスな技術が一番安全だし使えると思ってやっています。動かすときはLinuxでどうぞ。
コアダンプしないように作ってあるつもりです。（何らかの標準出力がある）
またブロックチェーン技術のためにネット上にあった自作でないSHA3を使い、電子署名モジュールは作るつもりです。
誰も作ってないようなので作りました。

実装は汚いコードですが、ポインタやオブジェクト指向を知らないので素人レベルのコードです。プロなら読めると思います。
全部で2000行くらいしかないです。
いつか勇者が現れて、私のコードを凌ぐ素晴らしいプログラムを作ってくれることを祈ります。
組み込み開発や、JavaScriptなんかにも直接移植できるのではないでしょうか。処理速度が遅いから無理かもしれないけど、
最適化は今後の課題です。まだプロセッサの能力をフルに発揮できていません。

残念ながらまだOpenMPには対応してません。
やり方が今一つよくわからないのですが、行列を並列計算させようとしても変な値になってしまいうまく並列化できません。
ちなみに開発環境は、LinuxでCPUはAMD　Ryzen7　2700Xです。マルチコアの意味がまるでないｗ。
ただその気になればパソコン用のCPUであれば、アセンブラマクロかなんかを使って、コア内で並列処理できると思います。
AMDのSIMDは遅いという評判ですがｗ。

参考文献を開発日記の別ファイルに移しました。今後ここにはマニュアル的なコードを書きたいと思います。

20200306（ひとことメモ）

機能拡張と計算結果の信頼性だけに注目しているので、全くコードを読む側のことを考えてません。
秘密鍵暗号を使って暗号化するのがいいか、ハッシュ関数を再帰的に使ってストリーム暗号のようにして使うかの比較。
実装が悪いのかあまりにも遅いので秘密鍵暗号やファイル暗号化はしないことに。

20200307

書くようなことでもないのだが、今のままでは鍵生成は遅く、並列化不可能。
AMD専用の最適化コンパイラclangを使って並列化できるところを並列化した。
自分の環境では33秒を切るくらいの時間で終わる。
特に用途はないので、ベンチマークにしてみてくださいｗ。
並列化オプションを外すと33秒後半。gccだと37秒位なので、いかにメーカー製のコンパイラが最適化しているか解る。
gccで速度が出ないということは基本的な設計から見直さないといけないのかもしれない。

20200318

いま高速化と最適化をしています。元からビックインディアンをリトルインディアンに直しているせいで鍵生成に無駄に時間がかかっています。
これは生成多項式が直感的に見やすくなるように次数の高い順から（配列の添字の小さい方から）代入をしているからで、これを改善すれば４倍くらい早くなりそうです。
といっても素人の実装なので期待しないでください。

20200320

最適化に挑戦してみました。バグトラップを外しでmakeで２７秒台を記録。


諸事情によりもう更新しません。
