# GoppaDecorder
耐量子計算機暗号のコアは、古典的Goppa符号の復号法にかかっているので、ここではその実装に役立つものを作ろうと思っています。
現在のところ、2の拡大体に関する復号は例外処理以外の実装で、バグがあります。
このコードを実装するのに使った書籍は、Oliver PrezzelのError-Correcting Codes and Finite Fieldsです。
実際に暗号を作るとなると、バイナリバージョンが必要なのですが、まだアルゴリズムを理解してないのでパターソンアルゴリズムを勉強中です。
このコードの最終目的は軍事レベルの暗号を作ることで、Niederreiter暗号を目指しています。鍵はでかいです。
毎月たくさんの新型が生まれている中で、なぜ改良版を使わないのかという人もいるかもしれませんが、
一番古いオーソドックスな技術が一番安全だし使えると思ってやっています。動かすときはLinuxでどうぞ。
コアダンプしないように作ってあるつもりです。（何らかの標準出力がある）
またブロックチェーン技術のためにネット上にあった自作でないSHA3を使い、電子署名モジュールは作るつもりです。
暗号だったらopensslがあるからいらないという人もいるかもしれないですが、Niederreiterは入ってません。GitHubも探しましたが誰も作ってないようです。
だったら私が、と重い腰をあげました。

実装は汚いコードですが、ポインタやオブジェクト指向を知らないので初心者でも読めるように書いてあるつもりです。
JavaScriptなんかにも直接移植できるのではないでしょうか。

メールはこちらまで：fumumue@gmail.com
