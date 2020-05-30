# GoppaDecorder(Goppie)

20200517

更新してないのに見に来てくれてありがとう。

数学をやるので、暫くの間更新しません。ブロックチェーンの実装をするためです。
では、またいつか。


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

20200321

特に変わってませんが、おすすめのコンパイラオプションを見つけたので書いておきます。
オリジナルのMakefileは汎用的なものですが、AMDのZENアーキテクチャに特化したバイナリコードを生成するためのオプションは次のようになります。

gcc -O3 -mtune=znver2 -march=znver2 -ffast-math -funroll-loops  -fopenmp oplib.c

これでバグトラップ有りで２６秒台がコンスタントに出ました。
gccもオプション次第でAMD純正コンパイラと遜色ない性能を出せることがわかりました。（ver 9.2）

最適化オプションを全くつけないと、１分２５秒でした。
素人の分際で最適化なんておこがましいと言われそうですが、何もしないよりマシだと思ったので、そのへんはコンパイラ任せにしようという感じですｗ．

20200322

使わない変数を少し削除した。
公開鍵計算のループを並列化に成功。なんと10秒で鍵生成から暗号化復号化までできる。
一方で、処理が高速化した関係か新たなバグが見つかった。
あまり喜んでもいられない。
よく考えてみれば鍵生成は最初の一回だけなので、20秒位は我慢できるかもしれない。
途中の処理でエラーがでるくらいなら、そのくらいの欠点には目をつぶるしかないのかもしれない。
並列化して生成された鍵を使うと5分もしないうちにエラーで止まってしまう。
こんなに頻繁にエラーが出ては使い物にならない。
デバッグしてもいいのだが、果たしてその努力に見合う価値があるかどうかは分からない。

20200223

並列処理で、鍵生成という一回限りの処理にバグがある。
でもそのバグを取る労力に見合う効果がないので、鍵生成の並列化はしないことにする。
今日は調子が悪いみたいでプログラムをしているとイライラする。
毎回鍵を生成しないようにファイルから読み込もうとしたのだがうまく行かないので、一昨日のプログラムに戻した。
そして一昨日のプログラムにバグが合ったのでそれを修正した。
差分の読み方がわからないから、ファイルビューで見てまるごとコピペしただけの一日。

追記：
言われるがままに分割コンパイルなんてなれないことをしたおかげで、昔動いていたファイルから鍵を読み込む方法がうまく行かなくなってしまった。
でも古いファイルをフォルダにまとめておいたのが動いたので事なきを得た。
我が道を行ったほうがいいようだ。
鍵生成なしだと暗号化復号化にOpenMPしないで9秒かかる。

20200324

バグを特定した。なんといつも決まった列が0になるというバグ。いつも418の倍数の列だけが0になる。
明らかな規則性を持つので、どこでそれが発生しているのかを調べればいいことになる。
これが明らかになれば、並列化の実現になることは間違いない。

原因判明。マルチスレッドの副作用。符号長６６８８を１６個に分割するとちょうど４１８になる。
１６個のスレッドが同時に動く時、一つのスレッドが４１８列を処理するのだが、ちょうど４１８番目の列が処理されずに残ってしまう。
つまり計算されない列が出てきてしまう。これを再計算して埋め合わさない限り完全な符号にはならない。
この手間を考えるとセキュリティ上完全な符号を作ったほうがいいので、わざわざ並列化して特定の列にエラーが来ないように小細工をするより馬鹿正直に１列づつ逐次処理をしたほうがいいということになる。
速さを取るか、セキュリティを取るかで考えれば２０秒の鍵生成のほうが確実だと言える。
短縮符号としてある特定の列を取り除き、計算しないという方法も考えられるが、そんなことをする意味がない。

追記：原因不明のバグもある。仮に４１８の倍数にエラーが来ないようにしても異常終了する場合がある。
その点、逐次処理をしている方には全くエラーが出ない。並列処理をしたものは５分と立たない間にエラーで止まる。

20200325

並列化しようとしても、変数が１つなので他のスレッドから上書きされて並列化出来ない。
最終的にほしいのは行列なので無理やり並列化する方法もあるかもしれないが、まだよく並列化アルゴリズムの動作原理を理解していない。
ループの処理したい範囲をスレッドに指定してやって、各スレッドごとに変数を用意してやればできる気がする。

追記；
para.cで生成多項式を固定した場合、パラレル処理とシリアル処理とで同じ鍵を生成することを確認した。
つまり同じ鍵をより早く並列処理で生成できたことになる。
しかしまだ油断できない。

20200326

並列化プログラミングでいくつかのキーワードを試したがうまく行かない。
全く同じ鍵ができるとは限らない。同じ鍵ができるときもあるという感じだ。スレッドIDがランダムに決まるせいかもしれない。
でももし違うとしてもシーケンシャルで生成した行列の列を入れ替えたものができるだけなので、本質的に同じ生成行列を出力していることになる。
そして実際そうやって作った生成行列も、もし同じならシーケンシャルと同じくらいの確実さで暗号化に使用できるはずだ。
というか、同じ鍵から生成したのではない生成行列だったら、そもそも復号が出来ないはずだから同じとみなしてもいいはずだ。

追記：
並列処理のテストバージョンである　test_prll.c　を追加。
ファイル内のdeta関数が行列の並列処理をして、高い確率でdetで逐次処理して生成した行列と一致する。
成功確率が高い気がするが、たまにハズレる。

20200327

引き続きtest_prll.cで耐久性実験を継続中。24時間動かしても止まらない。計算はかなり安定してきていると思う。
最も解空間が3000ビット以上もあるので1日動かしたからといって完璧だとは言えないのだが。
AMD Ryzen7 2700Xを25000円で買っておいてよかったと思う。
今のままでも十分早いし、買い替えるだけのメリットがないという評価だったので、Ryze3000シリーズはあまり評判のいいものではないようだ。
今買い替えるより次のプロセッサまで待った方がよさそう。去年買ったばかりなので減価償却してない。
時期プロセッサは一つのコアで処理できるスレッドの数が2倍になるという噂だ。

負荷のかかった状態で並列計算させた場合の結果が、AMD純正コンパイラとgccの間で差が出るようだ。
AMDのほうが早い気がするのは偶然だろうか。
例えば、グローバル変数の行列に複数のスレッドがアクセスしようとしたときに排他制御でブロックされたスレッドは書き込みに失敗して
その列が０になってしまうようだけど、AMDのコンパイラの実行ファイルだとこのような衝突を防いでくれるようでいつも一回でうまく行く。

20200328

引き続きでバック作業。確実視されていたdecode関数でエラーがでるので、pattarsonをコメントアウトして集中的にデバック作業を行っている。
一方で以下のサイトを見つけた。

http://arbitrage.jpn.org/it/ryzen-7-2700x/

このようにRyzenプロセッサはマルチスレッドが有利に作られているので、AMD純正コンパイラとマルチスレッドプログラミングは相性がいいと言っていい。
つまりOpenMPを多用するようなプログラミングで真の性能を引き出せるというわけだ。
Ryzenの性能を引き出すために、今後の目標としては量子計算機よりも並列計算の動作原理を中心に勉強するほうが大事らしい。

20200329

xgcdにバグを発見。とりあえず別関数を用意して、そちらで対処。decodeにはogcdを、pattartsonにはxgcdを使うことにしてデバッグ中。
一応終了条件を変えてみたけど、正しいかどうかは実験中。これでもまだエラーがでるならogcdを使うほうが安全かも。

20200330

24時間以上動いてますがまだエラーが出ません。やはりXGCDのバグだったのでしょうか？
並列計算で同じ鍵ができるというのは20連発成功という確率です。
今までNiederreiterの読み方がわからなかったんですが、ニーダーライダーというらしいです。かっこいいですね。

20200331

またエラーで停止。朝起きたら止まってた。気持ちが悪いことにxgcdの終了条件が２つになってしまう。
何か根本的に間違っているのだが、今はまだわからない。
decodeとpattarsonは終了条件が違うようなので、使用する関数をxgcdとogcdに分けました。

20200403

運用三日目。依然としてエラーは検出されない。これで完成と言ったところだろうか？
これからは電子署名と合わせたブロックチェーンを使ったスマートコントラクトを１から作ることになる。
PKI基盤や、ネット上に分散させるような秘密鍵の安全な管理に使う方向で行きたい。

20200404

実験終了。跡はバイナリ符号に置き換えるだけ。q進Goppaは完成だ。

20200405

https://onoono-life.blogspot.com/2014/04/ccopenmp.html

鍵生成の並列化にはこのサイトの情報が役に立ちました。
並列化に当たってはAMD純正のコンパイラを使いました。
gccで並列化を実行すると先に終了したスレッドがフライングをして変数のアクセスに失敗して何回もやり直す現象がおきます。
Makefileはgccでコンパイルできるように調整したものです。（まだ間違っているかもしれないので修正してください）
今日はバグのないバージョンにoplib.cを置き換えて、それ用にMakefileを書き換えました。
後はWallオプションで警告を消します。
バイナリ公開鍵を使えるようにします。

それにしてもコロナは一体いつ終息するんでしょうか。
これからもっとひどくなるんでしょうか・・・？

20200412

Makefileをgccで統一しました。
gccでも鍵を並列生成が可能です。
実験環境では約一分でした。

20200413

鍵生成をdetaにした時、正しく生成できないバグを発見。
現在デバッグ中。
clangでコンパイルしたときはエラー訂正に失敗しないのですが、gccでコンパイルすると１０分もしないでエラーで止まります。
しかも再現性のないエラーなので、gccのバグのような感じです。
ただ、aoccのclangではO3オプションが未実装なため、公開鍵を作る処理が遅く、最適化しても３分以上かかります。
これがgccだと１分で終わります。
そういうわけで完全な答えはないのですが、エラーが怖い人は鍵生成にdetaを使わないで、detを使うようにしてください。
鍵生成は別ファイルにして、aoccのclangでコンパイルした専用のコマンドにするかもしれないです。
大本のclangをコンパイルインストールしてみましたが、鍵生成は遅いままです。どうやらgcc -O3に当たるオプションは未実装のようです。
もし鍵生成を早くしたいのであればgccを指定してください。

20200416

clangによる鍵計算を高速化。
今まで3分30秒だった鍵計算処理を、pragma指定で2分短縮して1分40秒程度で完了できるようになった。
ちなみにgccは1分位。
そのまま実行するとスタックオーバーフローを起こすので、ulimit -s 16000と設定すること。
また、コンパイラはclangを使うことを強くおすすめする。
というのも、行列に複数のスレッドがアクセスしようとするときにclangは自動で排他制御をしているようなのだ。
一方gccの方は答えが一発で出ない。
繰り返して計算しているうちにclangと同じくらいの遅さになってしまう。
ここは正確さと速度の保証があるclangをおすすめしたい。
clangはaoccでなくても良いが、コンパイルにものすごく時間がかかるのでストレスになるかも。
でも計算結果は正確なのでおすすめする。

20200417

gccで鍵生成が1分4秒、clangで1分12秒位まで早くなった。
更新ファイルはlu.cのみ。

20200417

鍵生成の速度は昨日と変わらず。無理して早くすると正しい結果にならない。
今日は鍵生成から復号までの全体をチューニングした。

20200420

ダイエット中と体調不良により、暫くお休みします。
並列処理をするより、まず64ビット環境に最適化する必要があると思う。
何かいい方法が見つかったらやります。


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