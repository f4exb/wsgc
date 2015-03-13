

# Introduction #

We need a precise timing for the message. Firstly to identify the first PRN in a symbol as the message symbol is the repetition of at least 4 times the same PRN for proper integration. Secondly we also need to know the order in which the symbols are sent to be able to apply a Forward Error Correction scheme.

The proposed synchronization mechanism works with a training sequence allowing to make both synchronizations by identifying a point in time (time in the sense of the sequence of the PRNs in the message) at PRN resolution. The chip synchronization is already given by the correlation of the pilot.

It is assumed that a pilot sequence is sent in the signal with a PRN number different from the  pilot of messages to notify the receiver that a training sequence is being sent. The receiver will first look for a correlation with this pilot then when synchronizing properly it can switch to message processing on time at the start of the message sequence.

To allow receivers to synchronize (almost) randomly it is expected that training sequences are interspersed with CQ or beacon messages. Then when the QSO starts the training sequences can be omitted as with good performance clocks the receiver can stay in sync with the transmitter well beyond the end of the QSO. This is the general idea and details will be given when designing the exact QSO protocol.

The fun thing is that this scheme alone is capable of synchronizing two clocks at a distance (with the difference in propagation time) fairly accurately. If clocks are synchronized by another means like GPS then the propagation distance can be measured. The speedier the code the more precise but the more bandwidth needed.

# Summary of the idea #

The goal is to do an accumulation of correlation magnitudes with a sequence of PRNs rather than a repetition of PRNs as done with the message sequence.

The basic idea is to send a sequence of PRNs in the symbol alphabet order without repetition of PRNs. This order of PRNs is arbitrary and needs only to be known between transmitter and receiver. For one PRN period (synchronized by the pilot correlation) the results of the correlation are placed in the same PRN order and in an array of contiguous cells. For the next period the results are added shifted by one cell to the left (assuming horizontal representation and ascending order from left to right) thus the correlation values for PRNs sent add up on the same cell. After some integration time a peak should appear in the place of the first PRN received. Since this PRN is uniquely identified in the sequence (by its order in the alphabet) the point in time (in PRN periods) of the start of the training sequence can be determined. Then the message should start at a definite number of PRN periods from the start of the training sequence.

The full set of Gold Codes need not to be sent we only need to send a definite sequence among them.

This mechanism was validated with the C++ prototype and appears to work in conditions worse that the limit conditions for proper message decoding with signals affected with both AWGN and fading. Thus it is robust enough to ensure the synchronization of valid messages. An analysis window of at least 8 PRNs is recommended.

## A picture is better than a thousand words ##
This is a graphical representation of how the synchronization mechanism works. Correlation magnitudes for each PRN position are represented as vertical bars with a height respective to the magnitude value. Sequence was sent with PRN #0 first however the receiver catches up at PRN #3. The important point is that it is able to identify this was PRN #3.

![https://wsgc.googlecode.com/git/img/SyncScheme.png](https://wsgc.googlecode.com/git/img/SyncScheme.png)