package deconvolve

import (
	"container/heap"
	"fmt"
	"strings"
	"testing"
)

// The deadline list must hand back the entry furthest to the left first, and break ties by barcode name.
func TestEvictHeapOrder(t *testing.T) {
	var h evictHeap
	for _, e := range []evictEntry{{50, "b"}, {10, "z"}, {30, "a"}, {10, "a"}, {50, "a"}, {20, "m"}} {
		heap.Push(&h, e)
	}
	var got []string
	for h.Len() > 0 {
		e := heap.Pop(&h).(evictEntry)
		got = append(got, fmt.Sprintf("%d%s", e.end, e.bx))
	}
	if want := "10a 10z 20m 30a 50a 50b"; strings.Join(got, " ") != want {
		t.Errorf("pop order = %v, want %s", got, want)
	}
}
