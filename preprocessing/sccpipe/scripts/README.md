## Scripts docu
英文说明连自己都看不明白了，还是老老实实写点中文的吧
### polyA.pl
#### Mode1: 从左向右搜索，只要包含给定长度的A,就把这一段A以及右边的(3' prime)的序列都移除
#### Mode2: 在Mode1的基础上，加上了一个错配，也就是说会有更多的序列被删除
#### Mode3: 从右往左搜索，给定长度以及大于给定长度的A都会被删除。


``` shell
# you need to install a perl module before use it 
# type "cpan" in terminal then type: install Switch

Useage: perl polyA.pl input output badoutput polyA-length Mode

```
