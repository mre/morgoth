package ymtn

import (
	"bufio"
	"fmt"
)


//Output a given graph output into graphiz .dot format
func OutputGraph(out *bufio.Writer, wdg []*Node) {

	//Write header
	out.Write([]byte("digraph G {\n"))
	
	//Write all edges
	for _, node := range wdg {
		for _, edge := range node.Edges {
			out.Write([]byte(edge.Cause.Name()))
			out.Write([]byte(" -> "))
			out.Write([]byte(edge.Effect.Name()))
			out.Write([]byte(fmt.Sprintf("[label=\"%d\"];\n", edge.Weight)))
		}
	}

	//Write footer
	out.Write([]byte("}\n"))
	out.Flush()
}
