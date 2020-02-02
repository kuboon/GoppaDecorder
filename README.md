# GoppaDecorder

１．プログラム実行に当たって(oplib.c)　20200131

まず最初にmain関数のkey2()のコメントを外して鍵を作ります。鍵ができるまで最適化オプリョンを使うと３分程度で出来上がります。
次に復号のデモですがkey2()のコメントをまたつけなおしてください。すると、鍵ファイルから鍵を読み込んでランダムな位置に与えられたエラーの値、
この場合バイナリを意識しているので１しかないですが、それが８５個表示されます。暗号化と復号化に大体５秒位かかります。
今後はコメントやさらなるコードの可読性、速度を上げるための最適化、と言ってもASICなんかないので組み込み用というよりはパソコン用ですが、
などを保守していきたいです。組み込みで実行する際はフォークして改造してください。特に許可は必要ありません。ただし無保証です。
マルチコアの組み込みプロセッサなんかないだろうし、SIMDも使えないとなるともっと基本的な鍵生成アルゴリズムを見直さないといけないです。
スマホのプロセッサならマルチコアなので最適化して動かせるかもしれないです。ARMですね。開発環境がLinuxでVisualStudioもあるのですが、
どこまでいろんなプロセッサに対応できるかわかりません。スケルトンはできているはずなので、あとは保守だけです。
まあ使おうなんて言う人はいないでしょうが。

２．作るに当たった経緯

耐量子計算機暗号のコアの一つとして、古典的Goppa符号の復号法にかかっているので、ここではその実装に役立つものを作ろうと思っています。
現在のところ、2の拡大体に関する復号は例外処理以外の実装で、バグがあります。
このコードを実装するのに使った書籍は、Oliver PrezelのError-Correcting Codes and Finite Fieldsです。
実際に暗号を作るとなると、バイナリバージョンが必要なのですが、まだアルゴリズムを理解してないのでパターソンアルゴリズムを勉強中です。
このコードの最終目的は軍事レベルの暗号を作ることで、Niederreiter暗号を目指しています。鍵はでかいです。
毎月たくさんの新型が生まれている中で、なぜ改良版を使わないのかという人もいるかもしれませんが、
一番古いオーソドックスな技術が一番安全だし使えると思ってやっています。動かすときはLinuxでどうぞ。
コアダンプしないように作ってあるつもりです。（何らかの標準出力がある）
またブロックチェーン技術のためにネット上にあった自作でないSHA3を使い、電子署名モジュールは作るつもりです。
暗号だったらopensslがあるからいらないという人もいるかもしれないですが、Niederreiterは入ってません。
GitHubも探しましたがファイル一つで単体で動作するプログラムは誰も作ってないようです。
だったら私が、と重い腰をあげました。

実装は汚いコードですが、ポインタやオブジェクト指向を知らないので初心者レベルのコードです。プロなら読めると思います。
全部で１０００行くらいしかないです。
JavaScriptなんかにも直接移植できるのではないでしょうか。

こういう暗号のいいところは巨大整数を使ってないので簡単に組み込みレベルのプログラムに使えるという事です。
組み込み開発用に使うとしても、ほぼそのまま動くでしょう。残念ながらまだOpenMPには対応してません。
やり方が今一つよくわからないのですが、行列を並列計算させようとしても変な値になってしまいうまく並列化できません。
例えば今は鍵生成に長さ4096の列を170個一行ずつ計算して14分もかかっているのですが、これがマルチコアで計算できれば早くなることは間違いありません。
ちなみに開発環境は、LinuxでCPUはAMD　Ryzen7　2700Xです。マルチコアの意味がまるでないｗ。
ただその気になればパソコン用のCPUであれば、アセンブラマクロかなんかを使って、コア内で並列処理できると思います。
AMDのSIMDは遅いという評判ですがｗ。

３．開発日記

20191218

GF(4096)上定義されたバイナリGoppa符号を生成できるようになった。この計算には14分もかかるので今後何か別の計算方法を考えないといけない。(oplib.cpp)
因みに、QC-MDPCなんで邪道だと思っているので、オリジナルの古典的Goppa符号を使っている。

同時に、LU分解法にヒントを得た可逆スクランブル行列を生成できるようになった。(lu.c)

公開鍵のパラメーターは、[n,k,t]=[4096,2040,85]である。

20191221

GF4096において、型をunsigned shortにして、バイナリ生成行列を出力させる鍵生成まで完成。(oplib.cpp)

20191226

開発メモ。誤りロケータだけは実装完了。バイナリバージョンの場合は誤りの位置を決定するだけで済むのでこれで完成。(chen.c)

既知のバグ：誤りの値を計算するときに、0になる場合がある。これは誤りロケータと、誤り値関数にまだ最小公倍数、つまり
1次式が残されている証拠で、この点を解決するためにユークリッドアルゴリズムに多少の改良を加える必要がある。この点については
来年に目標にしたい。(oplib.c)

20191227

誤りロケータもバグっていることが判明。しかし体力の限界。もう十日近く15時間やってるけど、そろそろ限界。
今のところ一番確実に動くのがchen.cなので、これをもとにバグ取りをしようと思う。でも今年はもうおしまい。
わずか１０００行のプログラムの挙動が理解できない私は低能だｗ

手で計算できないものは計算機を使って計算させてバグ取りをしないといけないし、デバック用のプログラムも書かないといけないのでまた来年。

20191229

バグを最後の1つに絞り込めた。現在デバッグモジュールを作成中。

20191230

ついにバグを突き止めた。生成行列の特定の列が０ベクトルだったという落とし穴がｗ。アルゴリズムは間違いなかった。
どこを探しても問題なかったから、まさかと思って調べてみたら０だったというｗ。しかも１１個もあった。
陽性かくにん、よかった。あとは０のない生成行列をどう作るかで、トレースチェックだけすればいいんじゃないかと思うのだが違うんだろうか？
これでもかなりマシになったほうなのだが、まだ重解を持つ可能性があるので完全にバグがないとは言えない。
とりあえずバイナリの場合は誤り位置だけ正確に計算できればいいので、これで完成しているはず。
あとは実際に公開鍵を使って暗号化復号化をするというIOの部分を作ることになる。あとソースにコメントもつけないといけない。(oplib.c)

20191231

とりあえずバグはなくなったはずなのですが、まだあるかもしれません。そしてこのコードをどういう形で使うのかはまだ決めてません。
PGPのような完璧なアプリを作るのはちょっと荷が重いし、Cの初心者レベルのコードなので興味のある人が最適化したり、もっとエレガントな
コードを書いてくれればいいのですが。とりあえず趣味でやった割にはよくできたと自分を褒めてあげたいですｗ。

20200131

https://csrc.nist.gov/projects/post-quantum-cryptography/round-2-submissions

アメリカのNISTで進んでいる次世代公開鍵暗号の選定、第2ラウンドまでにMcEliece暗号が残っているようです。今の所符号に基づく暗号と講師に基づく暗号が拮抗しているようです。私もアセンブリマクロを使って最適化してみようかな。頑張れMcEliece！

https://klevas.mif.vu.lt/~skersys/vsd/crypto_on_codes/goppamceliece.pdf

入れ知恵になるかもしれないですがこれも参考にしてみたいと思います。どうも外部ライブラリを使った実装ばかりで私みたいに一枚岩の自己完結プログラムがないようで。

20200202

https://en.wikipedia.org/wiki/Binary_Goppa_code#Decoding

これのおかげでやっとパターソンアルゴリズムを理解できた。シンドローム多項式の逆多項式はゴッパ多項式が既約のとき、ユークリッドアルゴリズムで
計算できる。そしてバイナリGoppa符号の場合、符号長は変わらず次元だけが縦に引き伸ばされるので、シンドローム多項式はバイナリ多項式である。
このとき逆多項式とその平方が計算できるので、解くべき誤り位置関数の次数はシンドローム多項式の半分になる。だから拡大体のときより
Goppa多項式の次数は半分で済む。でも面倒だから今のままでもいいかなｗ

http://lupus.is.kochi-u.ac.jp/shiota/misc/field/FiniteField.html

